// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

/** \class
 *  \brief
 */


// =================================================================
//                              Includes
// =================================================================
#include <cassert>
#include <list>
#include <algorithm>
#include <memory>

#include "GeneTreeNode.h"
#include "GeneMutation.h"
#include "GeneticUnit.h"
#include "PointMutation.h"
#include "SmallInsertion.h"
#include "SmallDeletion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"


using std::list;


namespace aevol {

//##############################################################################
//                                                                             #
//                          Class GeneTreeNode                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
// =================================================================
//                             Constructors
// =================================================================
int32_t GeneTreeNode::nextID_ = 0;

GeneTreeNode::GeneTreeNode(int32_t nodeCreationDate, Protein * protein)
{
  ID_ = GeneTreeNode::nextID_;
  GeneTreeNode::nextID_ ++;

  dna_creation_date_ = nodeCreationDate;
  node_creation_date_ = nodeCreationDate;
  gene_loss_date_  = -1;
  gene_loss_type_  = NOT_LOST_YET;

  protein_pointer_ = protein;
  shine_dal_position_ = protein->shine_dal_pos();
  strand_ = protein->strand();

  nb_promoters_ = protein->rna_list().size();
  //printf("%d promoters at positions", nb_promoters_); // debug
  promoter_positions_ = new int32_t[nb_promoters_];
  rna_pointers_ = new Rna *[nb_promoters_];
  int32_t i = 0;
  for (const auto& rna: protein->rna_list()) {
    rna_pointers_[i] = rna;
    promoter_positions_[i] = rna->promoter_pos();
    i++;
  }
  // printf(" \n "); //debug

  left_child_  = NULL;
  right_child_ = NULL;
  parent_node_ = NULL;

  cds_possibly_modified_ = false;
  proms_possibly_modified_ = false;
  gene_possibly_duplicated_ = false;
  cds_completely_deleted_ = false;
  putative_position_for_the_duplicate_ = -1;
}


// =================================================================
//                             Destructors
// =================================================================

GeneTreeNode::~GeneTreeNode()
{
  if (left_child_ != NULL) delete left_child_;
  if (right_child_ != NULL) delete right_child_;
  if (promoter_positions_ != NULL) delete [] promoter_positions_;
  if (rna_pointers_ != NULL) delete [] rna_pointers_;
  std::for_each(mutation_list.begin(), mutation_list.end(), std::default_delete<GeneMutation>());
}



// =================================================================
//                            Public Methods
// =================================================================

GeneTreeNode *GeneTreeNode::search_in_subtree_leaves(const Protein * protein)
{
  GeneTreeNode *result_left = NULL, *result_right = NULL;
  if ((left_child_ == NULL) && (right_child_ == NULL)) // I am a leaf
    {
      if (protein_pointer_ == protein) return this;
      else return NULL;
    }
  else // I am an internal node
    {
      if (left_child_ != NULL)  result_left = left_child_->search_in_subtree_leaves(protein);
      if (right_child_ != NULL) result_right = right_child_->search_in_subtree_leaves(protein);
      if ((result_left == NULL) && (result_right == NULL)) return NULL;
      else if ((result_left == NULL) && (result_right != NULL)) return result_right;
      else if ((result_left != NULL) && (result_right == NULL)) return result_left;
      else
        {
          fprintf(stderr, "Error, the protein %p should not be found twice in the tree.\n", protein);
          // get the root_ address to print the whole tree to screen
          GeneTreeNode * n = parent_node_, *root = this;
          while (n!=NULL) {root = n; n = n->parent_node_; }
          // here, n==NULL and root points on the root of the tree
          root->print_subtree_to_screen();
          exit(EXIT_FAILURE);
        }
    }
}


void GeneTreeNode::print_subtree_to_screen()
{
  // Postorder tree traversal

  /* Left subtree */
  if (left_child_ != NULL)  left_child_->print_subtree_to_screen();

  /* Right subtree */
  if (right_child_ != NULL) right_child_->print_subtree_to_screen();

  /* Current tree node */
  printf("Node ID: %" PRId32 "\n", ID_);
  if (parent_node_ != NULL)   printf("Parent ID: %" PRId32 "\n", parent_node_->ID_);
  else                        printf("Parent ID: none\n");
  if (left_child_ != NULL)    printf("Left child ID: %" PRId32 "\n", left_child_->ID_);
  else                        printf("Left child ID: none\n");
  if (right_child_ != NULL)   printf("Right child ID: %" PRId32 "\n", right_child_->ID_);
  else                        printf("Right child ID: none\n");
  printf("Node creation date: %" PRId32 "\n", node_creation_date_);
  printf("DNA creation date: %" PRId32 "\n", dna_creation_date_);
  switch(gene_loss_type_)
    {
    case NOT_LOST_YET:            printf("Node status: NOT_LOST_YET\n"); break;
    case LOST_BY_LOCAL_MUTATION:  printf("Gene loss type: LOST_BY_LOCAL_MUTATION\n");  printf("Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case DELETED :                printf("Gene loss type: DELETED\n");  printf("Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case BROKEN_BY_REAR:          printf("Gene loss type: BROKEN_BY_REAR\n");  printf("Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case DUPLICATED:              printf("Node status: DUPLICATED\n");  printf("Duplication date: %" PRId32 "\n", gene_loss_date_);break;
    default: break;
    }
  printf("Protein pointer: %p, Shine-Dalgarno position: %" PRId32 "\n", protein_pointer_, shine_dal_position_);
  if(strand_ == LEADING) printf("Strand: LEADING\n");
  else                    printf("Strand: LAGGING\n");
  for(size_t i = 0; i < nb_promoters_; i++)
    {
      printf("Promoter at %" PRId32 ", rna pointer %p\n", promoter_positions_[i], rna_pointers_[i]);
    }
  printf("Number of mutations: %zu\n", mutation_list.size());
  for (const auto& mutation: mutation_list) {
    // TODO vld: simplify
    char str[128];
    mutation->description_string_for_gene_mut(str);
    printf("  %s\n", str);
  }
  printf("\n\n");
}


void GeneTreeNode::write_subtree_to_files(FILE * topologyFile, FILE * nodeAttributesFile, int32_t end_gener)
{
  // Newick format for the topology file (postorder tree traversal with parentheses and branch lengths)


  if ((left_child_ != NULL) || (right_child_ != NULL))
    {
      fprintf(topologyFile, "(");
      /* Left subtree */
      if (left_child_ != NULL) left_child_->write_subtree_to_files(topologyFile, nodeAttributesFile, end_gener);
     fprintf(topologyFile, ", ");
     /* Right subtree */
     if (right_child_ != NULL) right_child_->write_subtree_to_files(topologyFile, nodeAttributesFile, end_gener);
     fprintf(topologyFile, ")");
    }

  /* Current tree node */
  fprintf(topologyFile, "%" PRId32 "", ID_);
  if (gene_loss_type_ == NOT_LOST_YET)  fprintf(topologyFile, ":%" PRId32 "", end_gener - node_creation_date_);
  else fprintf(topologyFile, ":%" PRId32 "", gene_loss_date_ - node_creation_date_);


  fprintf(nodeAttributesFile, "Node ID: %" PRId32 "\n", ID_);
  if (parent_node_ != NULL)   fprintf(nodeAttributesFile, "Parent ID: %" PRId32 "\n", parent_node_->ID_);
  else                        fprintf(nodeAttributesFile, "Parent ID: none\n");
  if (left_child_ != NULL)    fprintf(nodeAttributesFile, "Left child ID: %" PRId32 "\n", left_child_->ID_);
  else                        fprintf(nodeAttributesFile, "Left child ID: none\n");
  if (right_child_ != NULL)   fprintf(nodeAttributesFile, "Right child ID: %" PRId32 "\n", right_child_->ID_);
  else                        fprintf(nodeAttributesFile, "Right child ID: none\n");
  fprintf(nodeAttributesFile, "Node creation date: %" PRId32 "\n", node_creation_date_);
  fprintf(nodeAttributesFile, "DNA creation date: %" PRId32 "\n", dna_creation_date_);
  switch(gene_loss_type_)
    {
    case NOT_LOST_YET:            fprintf(nodeAttributesFile, "Node status: NOT_LOST_YET\n"); break;
    case LOST_BY_LOCAL_MUTATION:  fprintf(nodeAttributesFile, "Gene loss type: LOST_BY_LOCAL_MUTATION\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case DELETED :                fprintf(nodeAttributesFile, "Gene loss type: DELETED\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case BROKEN_BY_REAR:          fprintf(nodeAttributesFile, "Gene loss type: BROKEN_BY_REAR\n"); fprintf(nodeAttributesFile, "Gene loss date: %" PRId32 "\n", gene_loss_date_); break;
    case DUPLICATED:              fprintf(nodeAttributesFile, "Node status:    DUPLICATED\n"); fprintf(nodeAttributesFile, "Duplication date: %" PRId32 "\n", gene_loss_date_); break;
    default: break;
    }
  if(strand_ == LEADING) fprintf(nodeAttributesFile, "Strand: LEADING\n");
  else                    fprintf(nodeAttributesFile, "Strand: LAGGING\n");
  fprintf(nodeAttributesFile, "Shine-Dalgarno position: %" PRId32 "\n", shine_dal_position_);
  for (size_t i= 0; i < nb_promoters_; i++)
    {
      fprintf(nodeAttributesFile, "Position of promoter %" PRId32 ": %" PRId32 "\n", static_cast<int32_t>(i+1), promoter_positions_[i]);
    }
  fprintf(nodeAttributesFile, "Number of mutations: %" PRId32 "\n", static_cast<int32_t>(mutation_list.size()));
  for (const auto& mutation: mutation_list) {
    // TODO vld: simplify
    char str[128];
    mutation->description_string_for_gene_mut(str);
    fprintf(nodeAttributesFile, "  %s", str);
    fprintf(nodeAttributesFile, "\n");
  }

  if (gene_loss_type_ == NOT_LOST_YET)
    {
      assert(protein_pointer_ != NULL);
      fprintf(nodeAttributesFile, "  Shine-Dalgarno pos:%" PRId32 ", Stop pos: %" PRId32 ", M: %.8f, W: %.8f, H: %.8f, nb promoters: %" PRId32 ", conc: %.8f \n", \
              protein_pointer_->shine_dal_pos(), protein_pointer_->last_STOP_base_pos(), \
              protein_pointer_->mean(), protein_pointer_->width(), protein_pointer_->height(),  \
              static_cast<int32_t>(protein_pointer_->rna_list().size()), protein_pointer_->concentration());
    }
  fprintf(nodeAttributesFile, "\n\n");
}



// all attributes on a single line
void GeneTreeNode::write_subtree_nodes_in_tabular_file(int32_t treeID, FILE * f)
{
  /* Left subtree */
  if (left_child_ != NULL) left_child_->write_subtree_nodes_in_tabular_file(treeID, f);
  /* Right subtree */
  if (right_child_ != NULL)  right_child_->write_subtree_nodes_in_tabular_file(treeID, f);


  fprintf(f, "%" PRId32 " ", treeID);
  fprintf(f, "%" PRId32 " ", ID_);
  if (parent_node_ != NULL)   fprintf(f, "%" PRId32 " ", parent_node_->ID_);
  else                        fprintf(f, "-1 ");
  if (left_child_ != NULL)    fprintf(f, "%" PRId32 " ", left_child_->ID_);
  else                        fprintf(f, "-1 ");
  if (right_child_ != NULL)   fprintf(f, "%" PRId32 " ", right_child_->ID_);
  else                        fprintf(f, "-1 ");
  fprintf(f, "%" PRId32 " ", node_creation_date_);
  fprintf(f, "%" PRId32 " ", dna_creation_date_);
  switch(gene_loss_type_)
    {
    case NOT_LOST_YET:            fprintf(f, "NOT_LOST_YET "); break;
    case LOST_BY_LOCAL_MUTATION:  fprintf(f, "LOST_BY_LOCAL_MUTATION "); break;
    case DELETED :                fprintf(f, "DELETED "); break;
    case BROKEN_BY_REAR:          fprintf(f, "BROKEN_BY_REAR "); break;
    case DUPLICATED:              fprintf(f, "DUPLICATED ");  break;
    default: break;
    }
  fprintf(f, "%" PRId32 " ", gene_loss_date_);
  if (strand_ == LEADING) fprintf(f, "LEADING ");
  else                    fprintf(f, "LAGGING ");
  fprintf(f, "%" PRId32 " ", shine_dal_position_);
  fprintf(f, "%" PRId32 " ", static_cast<int32_t>(nb_promoters_));
  if (gene_loss_type_ == NOT_LOST_YET)
    {
      assert(protein_pointer_ != NULL);
      fprintf(f, "%.8f %.8f %.8f %.8f ", \
              protein_pointer_->mean(), protein_pointer_->width(), protein_pointer_->height(),  \
              protein_pointer_->concentration());
    }
  else {fprintf(f, "-1 -1 -1 -1 ");}

  fprintf(f, "%" PRId32 " ", static_cast<int32_t>(mutation_list.size()));
  int32_t nb_localmut_upstream_neutral = 0, nb_localmut_upstream_benef = 0, nb_localmut_upstream_delet = 0;
  int32_t nb_rear_upstream_neutral = 0, nb_rear_upstream_benef = 0, nb_rear_upstream_delet = 0;
  int32_t nb_localmut_cds_neutral = 0, nb_localmut_cds_benef = 0, nb_localmut_cds_delet = 0;
  int32_t nb_rear_cds_neutral = 0, nb_rear_cds_benef = 0, nb_rear_cds_delet = 0;
  for (const auto& mutation: mutation_list) {
    if (mutation == &(*mutation_list.back()))
      // do not count the last event, if it was disruptive
      if ((gene_loss_type_ == DELETED) || (gene_loss_type_ == LOST_BY_LOCAL_MUTATION) || (gene_loss_type_ == BROKEN_BY_REAR))
        break;

    if (mutation->type_of_event() == 0) {
      if (mutation->region() == UPSTREAM) {
        if      (mutation->impact_on_metabolic_error() == 0.0) nb_localmut_upstream_neutral++;
        else if (mutation->impact_on_metabolic_error() < 0.0)  nb_localmut_upstream_benef++;
        else if (mutation->impact_on_metabolic_error() > 0.0)  nb_localmut_upstream_delet++;
      }
      else {
        if      (mutation->impact_on_metabolic_error() == 0.0) nb_localmut_cds_neutral++;
        else if (mutation->impact_on_metabolic_error() < 0.0)  nb_localmut_cds_benef++;
        else if (mutation->impact_on_metabolic_error() > 0.0)  nb_localmut_cds_delet++;
      }
    }
    else {
      if (mutation->region() == UPSTREAM) {
        if      (mutation->impact_on_metabolic_error() == 0.0) nb_rear_upstream_neutral++;
        else if (mutation->impact_on_metabolic_error() < 0.0)  nb_rear_upstream_benef++;
        else if (mutation->impact_on_metabolic_error() > 0.0)  nb_rear_upstream_delet++;
      }
      else {
        if      (mutation->impact_on_metabolic_error() == 0.0) nb_rear_cds_neutral++;
        else if (mutation->impact_on_metabolic_error() < 0.0)  nb_rear_cds_benef++;
        else if (mutation->impact_on_metabolic_error() > 0.0)  nb_rear_cds_delet++;
      }
    }
  }
  fprintf(f, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " ", \
          nb_localmut_upstream_neutral, nb_localmut_upstream_benef, nb_localmut_upstream_delet,  \
          nb_localmut_cds_neutral, nb_localmut_cds_benef, nb_localmut_cds_delet, \
          nb_rear_upstream_neutral, nb_rear_upstream_benef, nb_rear_upstream_delet, \
          nb_rear_cds_neutral, nb_rear_cds_benef, nb_rear_cds_delet);



  fprintf(f, "\n");

}



// This is an auxiliary function for the method anticipate_mutation_effect_on_genes_in_subtree.
// The segment should go from 'first' to 'last' (included) in the clockwise sense.
// 'first' and 'last' should not be equal.
static bool breakpoint_inside_segment(int32_t pos_brkpt, int32_t first, int32_t last)
{
  if (first < last) // most frequent case
    {
      if((first <= pos_brkpt) && (pos_brkpt <= last)) {return true;}
      else {return false;}
    }
  else // special case where the segment overlaps ori
    {
      if((first <= pos_brkpt) || (pos_brkpt <= last)) {return true;}
      else {return false;}
    }
}


// This is an auxiliary function for the method anticipate_mutation_effect_on_genes_in_subtree.
// It return true if the subsegment [first, last] is totally included in the segment [pos1, pos2].
// The subsegment should go from 'first' to 'last' in the clockwise sense and the segment should
// go from 'pos1' to 'pos2' in the clockwise sense.
static bool subsegment_totally_in_segment(int32_t pos1, int32_t pos2, int32_t first, int32_t last)
{
  if ((first < last)  && (pos1 <= pos2))
    {
      if (((first >= pos1) && (first <= pos2)) && ((last >= pos1) && (last <= pos2))) {return true; }
      else {return false;}
    }
  else if ((first < last) && (pos1 > pos2))  // mut seg in 2 pieces but not the gene
    {
      if ((first >= pos1) || (last <= pos2))  // the gene is either completely in [pos1, genlen-1] or completely in [0, pos2]
        {
          return true;
        }
      else return false;
    }
  else if ((first > last) && (pos1 <= pos2))  // gene in two pieces but not mut seg, the gene cannot be totally included
    {
      return false;
    }
  else // both mut seg and the gene are in 2 pieces
    {
      if ((first >= pos1) && (last <= pos2)) {return true;}
      else {return false;}
    }
}


void GeneTreeNode::update_pointers_in_subtree_leaves(GeneticUnit * unit)
{
 if ((left_child_ != NULL) || (right_child_ != NULL)) // I am a internal node
    {
      if (left_child_ != NULL)   left_child_->update_pointers_in_subtree_leaves(unit);
      if (right_child_ != NULL)  right_child_->update_pointers_in_subtree_leaves(unit);
    }
  else // no child => I am a leaf => there is work to do for me !
  {
    if (gene_loss_type_ != NOT_LOST_YET)
      // inactive leaf
      return;

    // TODO vld: refactor DUPLICATED CODE (ref dc1)
    auto& pl = unit->protein_list(strand_); // shorthand
    auto protein =
        find_if(pl.begin(), pl.end(),
                [this](Protein & p)
                {return p.shine_dal_pos() == shine_dal_position_;});
    if (protein != pl.end()) {
      /* The strand and shine dal position are correct */
      /* Update the protein and rna pointers and positions */
      nb_promoters_ = protein->rna_list().size();
      if (promoter_positions_ != NULL) delete [] promoter_positions_;
      if (rna_pointers_ != NULL) delete [] rna_pointers_;
      promoter_positions_ = new int32_t[nb_promoters_];
      rna_pointers_ = new Rna *[nb_promoters_];
      size_t i = 0;
      for (const auto& rna: protein->rna_list()) {
        rna_pointers_[i] = rna;
        promoter_positions_[i] = rna->promoter_pos();
        i++;
      }
    }
    else {
      fprintf(stderr, "Error: cannot find a protein that should be there.\n");
      exit(EXIT_FAILURE);
    }
  }
}

void GeneTreeNode::anticipate_mutation_effect_on_genes_in_subtree_leaves(const Mutation * mut, int32_t lengthOfGeneticUnit)
{
  if ((left_child_ != NULL) || (right_child_ != NULL)) // I am a internal node
    {
      if (left_child_ != NULL)   left_child_->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
      if (right_child_ != NULL)  right_child_->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
    }
  else // no child => I am a leaf => there is work to do for me !
    {
      if (gene_loss_type_ != NOT_LOST_YET)
        {
          // inactive leaf
          return;
        }

      int32_t genlen = lengthOfGeneticUnit; // in bp
      int32_t pos0 = -1, pos1 = -1, pos2 = -1, pos2bis = -1, pos3 = -1, mutlength = -1;
      // int32_t pos1donor = -1, pos2donor = -1, pos3donor = -1;  AlignmentSense sense = DIRECT;  // related to transfer (TO DO)
      bool invert = false;
      MutationType type = mut->mut_type();
      switch(type)
      {
        case SWITCH : {
          pos0 = dynamic_cast<const PointMutation*>(mut)->pos();
          mutlength = 1;
          break;
        }
        case S_INS : {
          const auto* s_ins = dynamic_cast<const SmallInsertion*>(mut);
          pos0 = s_ins->pos();
          mutlength = s_ins->length();
          break;
        }
        case S_DEL : {
          const auto* s_del = dynamic_cast<const SmallDeletion*>(mut);
          pos0 = s_del->pos();
          mutlength = s_del->length();
          break;
        }
        case DUPL : {
          const auto& dupl = dynamic_cast<const Duplication*>(mut);
          pos1 = dupl->pos1();
          pos2 = Utils::mod(dupl->pos2() - 1, genlen);
          pos2bis = dupl->pos2();
          pos0 = dupl->pos3();
          mutlength = dupl->length();
          break;
        }
        case DEL : {
          const auto& del = dynamic_cast<const Deletion*>(mut);
          pos1 = del->pos1();
          pos2 = Utils::mod(del->pos2() - 1, genlen);
          pos2bis = del->pos2();
          mutlength = del->length();
          break;
        }
        case TRANS : {
          const auto& trans = dynamic_cast<const Translocation*>(mut);
          pos1 = trans->pos1();
          pos2 = Utils::mod(trans->pos2() - 1, genlen);
          pos2bis = trans->pos2();
          pos3 = trans->pos3();
          pos0 = trans->pos4();
          invert = trans->invert();
          mutlength = trans->length();
          break;
        }
        case INV : {
          const auto& inv = dynamic_cast<const Inversion*>(mut);
          pos1 = inv->pos1();
          pos2 = Utils::mod(inv->pos2() - 1, genlen);
          pos2bis = inv->pos2();
          mutlength = inv->length();
          break;
        }
        case INSERT : {
          // TO DO
          break;
        }
        case INS_HT : {
          // TO DO
          break;
        }
        case REPL_HT : {
          // TO DO
          break;
        }
        default : {
          fprintf(stderr, "Error: unknown mutation type in GeneTreeNode::anticipate_mutation_effect_on_genes_in_subtree.\n");
        }
      }


      int32_t first_cds, last_cds;
      int32_t first_upstream, last_upstream; // "upstream region" is the segment between the furthest promoter and the Shine-Dalgarno sequence
      int32_t nbprom = protein_pointer_->rna_list().size();
      assert(nbprom != 0);
      assert(nbprom == static_cast<int32_t>(nb_promoters_));
      int32_t position_furthest_prom = -1, currentprompos = -1;
      if (protein_pointer_->strand() == LEADING)
        {
          first_cds = protein_pointer_->shine_dal_pos();
          last_cds = protein_pointer_->last_STOP_base_pos();
          for (const auto& rna: protein_pointer_->rna_list()) {
              currentprompos = rna->promoter_pos();
              if (currentprompos > first_cds) currentprompos = currentprompos - genlen; // negative value for promoters on the other side of ori
              if ((position_furthest_prom == -1) || (position_furthest_prom < currentprompos)) // we need the smallest promoter position
                {
                  position_furthest_prom = rna->promoter_pos();
                }
          }
          position_furthest_prom = Utils::mod(position_furthest_prom, genlen); // restore a positive value if necessary
          first_upstream = position_furthest_prom;
          last_upstream = Utils::mod(first_cds - 1, genlen);
        }
      else
        {
          first_cds = protein_pointer_->last_STOP_base_pos();
          last_cds = protein_pointer_->shine_dal_pos();
          for (const auto& rna: protein_pointer_->rna_list()) {
              currentprompos = rna->promoter_pos();
              if (currentprompos < last_cds) currentprompos = currentprompos + genlen; // value larger than genlen for promoters on the other side of ori
              if ((position_furthest_prom == -1) || (position_furthest_prom > currentprompos)) // we need the largest promoter position
                {
                  position_furthest_prom = rna->promoter_pos();
                }
          }
          position_furthest_prom = Utils::mod(position_furthest_prom, genlen); // restore a value < genlen if necessary
          first_upstream = Utils::mod(last_cds - 1, genlen);
          last_upstream = position_furthest_prom;
        }




      switch(type)
        {
        case SWITCH:
          {
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) cds_possibly_modified_ = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) proms_possibly_modified_ = true;
            break;
          }
        case S_INS:
          {
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) cds_possibly_modified_ = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) proms_possibly_modified_ = true;

            if (shine_dal_position_ >= pos0) shine_dal_position_ = Utils::mod(shine_dal_position_ + mutlength, genlen + mutlength);
            for (size_t i = 0; i < nb_promoters_; i++)
              {
                if (promoter_positions_[i] >= pos0){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] + mutlength, genlen + mutlength);}
              }

            break;
          }
        case S_DEL:
          {
            // If the Shine-Dalgarno position is in the small deleted segment, mark the cds as deleted
            // (our gene tracking is based on the tracking of the Shine-Dalgarno position,
            // and we cannot predict the position of a bp that was deleted: we lose track of the gene)
            if (mutlength == 1)
              {
                if (protein_pointer_->shine_dal_pos() == pos0)  cds_completely_deleted_ = true;
              }
            else // mutlength > 1
              {
                if (breakpoint_inside_segment(protein_pointer_->shine_dal_pos(), pos0, Utils::mod(pos0 + mutlength - 1, genlen))) cds_completely_deleted_ = true;
              }

            if (!(cds_completely_deleted_))
              {
                if (breakpoint_inside_segment(pos0, Utils::mod(first_cds - mutlength, genlen), last_cds)) cds_possibly_modified_ = true;
                if (breakpoint_inside_segment(pos0, Utils::mod(first_upstream - mutlength, genlen), last_upstream)) proms_possibly_modified_ = true;

                if (pos0 + mutlength <= genlen) // the deletion does not contain the replication origin
                  {
                    if (shine_dal_position_ >= pos0) shine_dal_position_ = Utils::mod(shine_dal_position_ - mutlength, genlen - mutlength);
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if (promoter_positions_[i] >= pos0){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] - mutlength, genlen - mutlength);}
                      }
                  }
                else // the deletion contains the replication origin
                  {
                    int32_t nb_del_after_ori = mutlength - genlen + pos0;
                    if (shine_dal_position_ >= 0) shine_dal_position_ = Utils::mod(shine_dal_position_ - nb_del_after_ori, genlen - mutlength);
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if (promoter_positions_[i] >= 0){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] - nb_del_after_ori, genlen - mutlength);}
                      }
                  }
              }
            break;
          }
        case DUPL:
          {
            if (subsegment_totally_in_segment(pos1, pos2, first_cds, last_cds))
              {
                gene_possibly_duplicated_ = true;
                putative_position_for_the_duplicate_ = Utils::mod(
                    Utils::mod(shine_dal_position_ - pos1, genlen) + pos0, genlen + mutlength);
              }
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) cds_possibly_modified_ = true;
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) proms_possibly_modified_ = true;

            if (shine_dal_position_ >= pos0) shine_dal_position_ = Utils::mod(shine_dal_position_ + mutlength, genlen + mutlength);
            for (size_t i = 0; i < nb_promoters_; i++)
              {
                if (promoter_positions_[i] >= pos0){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] + mutlength, genlen + mutlength);}
              }

            break;
          }
        case DEL:
          {
            // If the Shine-Dalgarno is in the deleted segment, mark the cds as deleted
            // (our gene tracking is based on the tracking of the Shine-Dalgarno position,
            // and we cannot predict the position of a bp that was deleted: we lose track of the gene)
           if (mutlength == 1)
              {
                if (protein_pointer_->shine_dal_pos() == pos1)  cds_completely_deleted_ = true;
              }
            else // mutlength > 1
              {
                if (breakpoint_inside_segment(protein_pointer_->shine_dal_pos(), pos1, pos2)) cds_completely_deleted_ = true;
              }

            if (!(cds_completely_deleted_))
              {
                if (breakpoint_inside_segment(pos1, first_cds, last_cds) || breakpoint_inside_segment(pos2, first_cds, last_cds)) cds_possibly_modified_ = true;
                if (breakpoint_inside_segment(pos1, first_upstream, last_upstream) || breakpoint_inside_segment(pos2, first_upstream, last_upstream)) proms_possibly_modified_ = true;

                if (pos1 < pos2bis) // the deletion does not contain the replication origin
                  {
                    if (shine_dal_position_ >= pos1) shine_dal_position_ = Utils::mod(shine_dal_position_ - mutlength, genlen - mutlength);
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if (promoter_positions_[i] >= pos1){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] - mutlength, genlen - mutlength);}
                      }
                  }
                else  // the deletion contains the replication origin
                  {
                    if (shine_dal_position_ >= 0) shine_dal_position_ = Utils::mod(shine_dal_position_ - pos2bis, genlen - mutlength);
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if (promoter_positions_[i] >= 0){ promoter_positions_[i] = Utils::mod(promoter_positions_[i] - pos2bis, genlen - mutlength);}
                      }
                  }
              }
            break;
          }
        case TRANS:
          {
            if (breakpoint_inside_segment(pos1, first_cds, last_cds)) cds_possibly_modified_ = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_cds, last_cds)) cds_possibly_modified_ = true;  // end of the excised segment
            if (breakpoint_inside_segment(pos3, first_cds, last_cds)) cds_possibly_modified_ = true;  // breakpoint inside the segment for the reinsertion
            if (breakpoint_inside_segment(pos0, first_cds, last_cds)) cds_possibly_modified_ = true;  // reinsertion point in the genetic unit
            if (breakpoint_inside_segment(pos1, first_upstream, last_upstream)) proms_possibly_modified_ = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_upstream, last_upstream)) proms_possibly_modified_ = true;  // end of the excised segment
            if (breakpoint_inside_segment(pos3, first_upstream, last_upstream)) proms_possibly_modified_ = true;  // breakpoint inside the segment for the reinsertion
            if (breakpoint_inside_segment(pos0, first_upstream, last_upstream)) proms_possibly_modified_ = true;  // reinsertion point in the genetic unit

            int32_t pos_min = Utils::min(pos1, Utils::min(pos2bis, Utils::min(pos3, pos0)));
            int32_t pos_B, pos_C, pos_D, pos_E;
            int32_t len_B, len_C, len_D;
            if (! invert)
              {
                if (pos_min == pos1)          { pos_B = pos1; pos_C = pos3; pos_D = pos2bis; pos_E = pos0; }
                else if (pos_min == pos2bis)  { pos_B = pos2bis; pos_C = pos0; pos_D = pos1; pos_E = pos3; }
                else if (pos_min == pos3)     { pos_B = pos3; pos_C = pos2bis; pos_D = pos0; pos_E = pos1; }
                else                            { pos_B = pos0; pos_C = pos1; pos_D = pos3; pos_E = pos2bis; } // if (pos_min == pos0)
                len_B = pos_C - pos_B;
                len_C = pos_D - pos_C;
                len_D = pos_E - pos_D;
                if      ((shine_dal_position_ >= pos_B) && (shine_dal_position_ < pos_C))   shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D + len_C, genlen);
                else if ((shine_dal_position_ >= pos_C) && (shine_dal_position_ < pos_D))   shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D - len_B, genlen);
                else if ((shine_dal_position_ >= pos_D) && (shine_dal_position_ < pos_E))   shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B - len_C, genlen);
                for (size_t i = 0; i < nb_promoters_; i++)
                  {
                    if      ((promoter_positions_[i] >= pos_B) && (promoter_positions_[i] < pos_C))  promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D + len_C, genlen);
                    else if ((promoter_positions_[i] >= pos_C) && (promoter_positions_[i] < pos_D))  promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D - len_B, genlen);
                    else if ((promoter_positions_[i] >= pos_D) && (promoter_positions_[i] < pos_E))  promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B - len_C, genlen);
                  }
              }
            else // invert
              {
                if (pos_min == pos1)
                  {
                    pos_B = pos1; pos_C = pos3; pos_D = pos2bis; pos_E = pos0;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((shine_dal_position_ >= pos_B) && (shine_dal_position_ < pos_C))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_B + pos_C - shine_dal_position_ -1 ;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_C) && (shine_dal_position_ < pos_D))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_C + pos_D - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_D) && (shine_dal_position_ < pos_E))
                      {
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B - len_C, genlen);
                      }
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if  ((promoter_positions_[i] >= pos_B) && (promoter_positions_[i] < pos_C))
                          {
                            promoter_positions_[i] = pos_B + pos_C - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_C) && (promoter_positions_[i] < pos_D))
                          {
                            promoter_positions_[i] = pos_C + pos_D - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_D) && (promoter_positions_[i] < pos_E))
                          {
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B - len_C, genlen);
                          }
                      }
                  }
                else if (pos_min == pos2bis)
                  {
                    pos_B = pos2bis; pos_C = pos0; pos_D = pos1; pos_E = pos3;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((shine_dal_position_ >= pos_B) && (shine_dal_position_ < pos_C))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_B + pos_C - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_C) && (shine_dal_position_ < pos_D))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_C + pos_D - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_D) && (shine_dal_position_ < pos_E))
                      {
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B - len_C, genlen);
                      }
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if  ((promoter_positions_[i] >= pos_B) && (promoter_positions_[i] < pos_C))
                          {
                            promoter_positions_[i] = pos_B + pos_C - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_C) && (promoter_positions_[i] < pos_D))
                          {
                            promoter_positions_[i] = pos_C + pos_D - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_D) && (promoter_positions_[i] < pos_E))
                          {
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B - len_C, genlen);
                          }
                      }
                  }
                else if (pos_min == pos3)
                  {
                    pos_B = pos3; pos_C = pos2bis; pos_D = pos0; pos_E = pos1;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((shine_dal_position_ >= pos_B) && (shine_dal_position_ < pos_C))
                      {
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_C + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_C) && (shine_dal_position_ < pos_D))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_C + pos_D - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_D) && (shine_dal_position_ < pos_E))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_D + pos_E - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B, genlen);
                      }
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if  ((promoter_positions_[i] >= pos_B) && (promoter_positions_[i] < pos_C))
                          {
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_C + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_C) && (promoter_positions_[i] < pos_D))
                          {
                            promoter_positions_[i] = pos_C + pos_D - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_D) && (promoter_positions_[i] < pos_E))
                          {
                            promoter_positions_[i] = pos_D + pos_E - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B, genlen);
                          }
                      }

                  }
                else // if (pos_min == pos0)
                  {
                    pos_B = pos0; pos_C = pos1; pos_D = pos3; pos_E = pos2bis;
                    len_B = pos_C - pos_B;
                    len_C = pos_D - pos_C;
                    len_D = pos_E - pos_D;
                    if  ((shine_dal_position_ >= pos_B) && (shine_dal_position_ < pos_C))
                      {
                        shine_dal_position_ = Utils::mod(shine_dal_position_ + len_C + len_D, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_C) && (shine_dal_position_ < pos_D))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_C + pos_D - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B, genlen);
                      }
                    else if ((shine_dal_position_ >= pos_D) && (shine_dal_position_ < pos_E))
                      {
                        if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                        shine_dal_position_ = pos_D + pos_E - shine_dal_position_ - 1;
                        shine_dal_position_ = Utils::mod(shine_dal_position_ - len_B, genlen);
                      }
                    for (size_t i = 0; i < nb_promoters_; i++)
                      {
                        if  ((promoter_positions_[i] >= pos_B) && (promoter_positions_[i] < pos_C))
                          {
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] + len_C + len_D, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_C) && (promoter_positions_[i] < pos_D))
                          {
                            promoter_positions_[i] = pos_C + pos_D - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B, genlen);
                          }
                        else if ((promoter_positions_[i] >= pos_D) && (promoter_positions_[i] < pos_E))
                          {
                            promoter_positions_[i] = pos_D + pos_E - promoter_positions_[i] - 1;
                            promoter_positions_[i] = Utils::mod(promoter_positions_[i] - len_B, genlen);
                          }
                      }
                  }
              }

            break;
          }
        case INV:
          {
            if (breakpoint_inside_segment(pos1, first_cds, last_cds)) cds_possibly_modified_ = true;
            if (breakpoint_inside_segment(pos2, first_cds, last_cds)) cds_possibly_modified_ = true;
            if (breakpoint_inside_segment(pos1, first_upstream, last_upstream)) proms_possibly_modified_ = true;   // beginning of the excised segment
            if (breakpoint_inside_segment(pos2, first_upstream, last_upstream)) proms_possibly_modified_ = true;  // end of the excised segment

            if  ((shine_dal_position_ >= pos1) && (shine_dal_position_ < pos2bis))
              {
                if (strand_ == LEADING) {strand_ = LAGGING;} else {strand_ = LEADING;}
                shine_dal_position_ = Utils::mod(pos1 + pos2bis - shine_dal_position_ - 1, genlen);
              }
            break;
          }
        case INSERT:
          {
            // TO DO
            break;
          }
        case INS_HT:
          {
            // TO DO
            break;
          }
        case REPL_HT:
          {
            // TO DO
            break;
          }
        default:
          // Only simple mutation types are considered.
          break;
        }
    }
}




void GeneTreeNode::register_actual_mutation_effect_on_genes_in_subtree_leaves(
    GeneTree * tree, const Mutation * mut, GeneticUnit* unit, int32_t gener, double impact_on_metabolic_error)
{
  if ((left_child_ != NULL) || (right_child_ != NULL)) // I am a internal node, just delegate work to others
    {
      if (left_child_ != NULL)   left_child_->register_actual_mutation_effect_on_genes_in_subtree_leaves(tree, mut, unit, gener, impact_on_metabolic_error);
      if (right_child_ != NULL)  right_child_->register_actual_mutation_effect_on_genes_in_subtree_leaves(tree, mut, unit, gener, impact_on_metabolic_error);
    }
  else // no child => I am a leaf => there is work to do for me !
    {

      if (gene_loss_type_ != NOT_LOST_YET)
        {
          // inactive leaf
          return;
        }

      GeneMutation * genemut = NULL;

      if (cds_completely_deleted_)
        {
          genemut = new GeneMutation(*mut, gener, shine_dal_position_, strand_, CDS);
          genemut->set_impact_on_metabolic_error(impact_on_metabolic_error);
          mutation_list.push_back(genemut);
          if (gener > tree->end_gener_) tree->end_gener_ = gener;

          gene_loss_date_ = gener;
          gene_loss_type_ = DELETED;
          protein_pointer_ = NULL;
          for (size_t i = 0; i < nb_promoters_; i++) {rna_pointers_[i] = NULL;}
          if (gener > tree->end_gener_) tree->end_gener_ = gener;
          (tree->nb_active_leaves_) --;
          return;
        }


      if ((!cds_completely_deleted_) && (!cds_possibly_modified_) && (!proms_possibly_modified_))
        {
          // This CDS was not affected by the mutation (it could have be moved or duplicated however).
          // Just make sure that we have correctly predicted the positions of the SD sequence and of the promoters.

          // TODO vld: refactor DUPLICATED CODE (ref dc1)
          auto& pl = unit->protein_list(strand_);
          auto protein =
              find_if(pl.begin(), pl.end(),
                      [this](Protein & p)
                      { return p.shine_dal_pos() == shine_dal_position_; });
          if (protein != pl.end()) {
            /* The strand and shine dal position are correct */
            /* Update the protein and rna pointers and positions */
            nb_promoters_ = protein->rna_list().size();
            if (promoter_positions_ != NULL) delete [] promoter_positions_;
            if (rna_pointers_ != NULL) delete [] rna_pointers_;
            promoter_positions_ = new int32_t[nb_promoters_];
            rna_pointers_ = new Rna *[nb_promoters_];
            size_t i = 0;
            for (const auto& rna: protein->rna_list()) {
              rna_pointers_[i] = rna;
              promoter_positions_[i] = rna->promoter_pos();
              i++;
            }
          }
          else
            {
              fprintf(stderr, "Error: cannot find a protein that should have survived.\n");
              char str[100];
              mut->generic_description_string(str);
              printf("Mutation : %s\n", str);
              printf("CDS should be at %d ", shine_dal_position_);
              if (strand_ == LEADING) printf("LEADING\n");
              else printf("LAGGING\n");
              unit->print_proteins();
              exit(EXIT_FAILURE);
            }
        }


      if (cds_possibly_modified_ || proms_possibly_modified_)
        {
          /* Record the impact of the mutation on the metabolic error */
          if      ((cds_possibly_modified_) && (proms_possibly_modified_))   genemut = new GeneMutation(*mut, gener, shine_dal_position_, strand_, BOTH);
          else if ((cds_possibly_modified_) && (!proms_possibly_modified_))  genemut = new GeneMutation(*mut, gener, shine_dal_position_, strand_, CDS);
          else if ((!cds_possibly_modified_) && (proms_possibly_modified_))  genemut = new GeneMutation(*mut, gener, shine_dal_position_, strand_, UPSTREAM);
          genemut->set_impact_on_metabolic_error(impact_on_metabolic_error);
          mutation_list.push_back(genemut);
          if (gener > tree->end_gener_) tree->end_gener_ = gener;

          /* Check whether the protein survived the event */
          // TODO vld: refactor DUPLICATED CODE (ref dc1)
          auto& pl = unit->protein_list(strand_); // shorthand
          auto protein =
              find_if(pl.begin(), pl.end(),
                      [this](Protein & p)
                      { return p.shine_dal_pos() == shine_dal_position_; });
          if (protein != pl.end()) {
            /* The strand and shine dal position are correct */
            /* Update the protein and rna pointers and positions */
            nb_promoters_ = protein->rna_list().size();
            if (promoter_positions_ != NULL) delete [] promoter_positions_;
            if (rna_pointers_ != NULL) delete [] rna_pointers_;
            promoter_positions_ = new int32_t[nb_promoters_];
            rna_pointers_ = new Rna *[nb_promoters_];
            size_t i = 0;
            for (const auto& rna: protein->rna_list()) {
              rna_pointers_[i] = rna;
              promoter_positions_[i] = rna->promoter_pos();
              i++;
            }
          }
          else
            {
              /* The protein does not exist anymore, the gene was killed by the event */
              gene_loss_date_ = gener;
              if ((mut->mut_type() == SWITCH) || (mut->mut_type() == S_INS) || (mut->mut_type() == S_DEL)) gene_loss_type_ = LOST_BY_LOCAL_MUTATION;
              else if ((mut->mut_type() == DUPL) || (mut->mut_type() == DEL) || (mut->mut_type() == TRANS) || (mut->mut_type() == INV)) gene_loss_type_ = BROKEN_BY_REAR;
              protein_pointer_ = NULL;
              for (size_t i = 0; i < nb_promoters_; i++) {rna_pointers_[i] = NULL;}
              if (gener > tree->end_gener_) (tree->end_gener_) = gener;
              (tree->nb_active_leaves_) --;

            }

        }


      if (gene_possibly_duplicated_) {
        /* Check whether the duplicated CDS found a promoter */
        /* It should be on the same strand as myself, at the putative_position_for_the_duplicate_ */

        auto& pl = unit->protein_list(strand_); // shorthand
        auto protein =
            find_if(pl.begin(), pl.end(),
                    [this](Protein & p)
                    { return p.shine_dal_pos() == putative_position_for_the_duplicate_; });

        if (protein != pl.end()) {
          if (protein_pointer_ != NULL) {
            // Create a new node for the "old" DNA segment
            left_child_ = new GeneTreeNode(gener, protein_pointer_);
            left_child_->node_creation_date_ = gener;
            left_child_->dna_creation_date_ = dna_creation_date_;
            left_child_->parent_node_ = this;
          }
          else {
            // Poor old gene was killed by the insertion of the duplicated segment
            // We do not create a left child
            left_child_ = NULL;
          }

          // Create a new node for the "new" DNA segment
          right_child_ = new GeneTreeNode(gener, &*protein);
          right_child_->node_creation_date_ = gener;
          right_child_->dna_creation_date_ = gener;
          right_child_->parent_node_ = this;

          // This node becomes internal, it represents an ancestral (obsolete) state of the gene
          protein_pointer_ = NULL;
          for (size_t i = 0; i < nb_promoters_; i++)
            rna_pointers_[i] = NULL;
          gene_loss_type_ = DUPLICATED;
          gene_loss_date_ = gener;

          // Update tree statistics
          if (left_child_ != NULL) {
            (tree->total_nb_nodes_) += 2;
            (tree->nb_internal_nodes_) ++;
            (tree->nb_leaves_) ++;  // - 1 + 2 (the ex-leaf becomes an internal node, 2 leaves are created)

          }
          else {
            (tree->total_nb_nodes_) += 1;
            (tree->nb_internal_nodes_) ++;
            // (tree->nb_leaves_) remains unchanged  <==  - 1 + 1 (the ex-leaf becomes an internal node, 1 leave is created)
          }
          if (protein != pl.end())
            tree->nb_active_leaves_++;
          if (gener > tree->end_gener_)
            tree->end_gener_ = gener;
        }
        // else nothing to do, the duplication was only partial, not a complete gene duplication
      }


      /* Get ready for the next mutation */
      cds_possibly_modified_ = false;
      proms_possibly_modified_ = false;
      gene_possibly_duplicated_ = false;
      cds_completely_deleted_ = false;
      putative_position_for_the_duplicate_ = -1;
    }

}



// =================================================================
//                           Protected Methods
// =================================================================


} // namespace aevol
