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
//                              Libraries
// =================================================================
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>



// =================================================================
//                            Project Files
// =================================================================

#include "GeneTree.h"
#include "GeneTreeNode.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                            Class GeneTree                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================



// =================================================================
//                             Constructors
// =================================================================

GeneTree::GeneTree()
{
  root_ = NULL;
  begin_gener_ = 0;
  end_gener_ = 0;
  total_nb_nodes_ = 0;
  nb_internal_nodes_ = 0;
  nb_leaves_ = 0;
  nb_active_leaves_ = 0;
  creation_type_ = INITIALIZATION;
}

// Creates a tree with just a root node.
GeneTree::GeneTree(int32_t nodeCreationDate, Protein * protein, const Mutation * mut /* = NULL */)
{
  root_ = new GeneTreeNode(nodeCreationDate, protein);
  if (mut == NULL)  creation_type_ = INITIALIZATION;
  else if ((mut->mut_type() == SWITCH) || (mut->mut_type() == S_INS) || (mut->mut_type() == S_DEL)) creation_type_ = LOCAL_MUTATION;
  else if ((mut->mut_type() == DUPL) || (mut->mut_type() == DEL) || (mut->mut_type() == TRANS) || (mut->mut_type() == INV)) creation_type_ = REARRANGEMENT;
  else creation_type_ = TRANSFER;


  begin_gener_ = nodeCreationDate;
  end_gener_ = nodeCreationDate;
  total_nb_nodes_ = 1;
  nb_internal_nodes_ = 0;
  nb_leaves_ = 1;
  if (protein != NULL) nb_active_leaves_ = 1; else nb_active_leaves_ = 0;
}




// =================================================================
//                             Destructors
// =================================================================

GeneTree::~GeneTree()
{
  delete root_;
}


// =================================================================
//                           Public Methods
// =================================================================

void GeneTree::set_end_gener_if_active_leaves(int32_t gener)
{
  if (nb_active_leaves_ > 0) end_gener_ = gener;
}

void GeneTree::update_pointers_in_tree_leaves(GeneticUnit * unit)
{
  root_->update_pointers_in_subtree_leaves(unit);
}

void GeneTree::anticipate_mutation_effect_on_genes_in_tree_leaves(const Mutation * mut, int32_t lengthOfGeneticUnit)
{
  root_->anticipate_mutation_effect_on_genes_in_subtree_leaves(mut, lengthOfGeneticUnit);
}

void GeneTree::register_actual_mutation_effect_on_genes_in_tree_leaves(const Mutation * mut, GeneticUnit * unit, int32_t gener, double impact_on_metabolic_error)
{
  root_->register_actual_mutation_effect_on_genes_in_subtree_leaves(this, mut, unit, gener, impact_on_metabolic_error);
}

// void GeneTree::duplicate_this_gene(GeneTreeNode * node, int32_t duplicDate, Protein * newProtein)
// {
//   if (newProtein == node->protein_pointer_) {fprintf(stderr, "Error, duplication with the same protein\n"); exit(EXIT_FAILURE);}

//   // Create a new node for the "old" DNA segment
//   node->left_child_ = new GeneTreeNode(duplicDate, node->protein_pointer_);
//   node->left_child_->node_creation_date_ = duplicDate;
//   node->left_child_->dna_creation_date_ = node->dna_creation_date_;
//   node->left_child_->parent_node_ = node;

//   // Create a new node for the "new" DNA segment
//   node->right_child_ = new GeneTreeNode(duplicDate, newProtein);
//   node->right_child_->node_creation_date_ = duplicDate;
//   node->right_child_->dna_creation_date_ = duplicDate;
//   node->right_child_->parent_node_ = node;

//   // This node becomes internal, it represents an ancestral (obsolete) state of the gene
//   node->protein_pointer_ = NULL;
//   for (int32_t i = 0; i < node->nb_promoters_; i++) {node->rna_pointers_[i] = NULL;}
//   node->gene_loss_type_ = DUPLICATED;
//   node->gene_loss_date_ = duplicDate;

//   // Update tree statistics
//   total_nb_nodes_ += 2;
//   nb_internal_nodes_ ++;
//   nb_leaves_ ++;  // - 1 + 2 (the ex-leaf becomes an internal node, 2 leaves are created)
//   if (newProtein != NULL) nb_active_leaves_ ++;
//   if (duplicDate > end_gener_) end_gener_ = duplicDate;
// }


// void GeneTree::report_gene_mutation(GeneTreeNode * node, GeneMutation * geneMut)
// {
//   node->mutation_list_->add(geneMut);
//   if (geneMut->generation() > end_gener_) end_gener_ = geneMut->generation();
// }



// void GeneTree::report_gene_loss(GeneTreeNode * node, int32_t geneLossDate, ae_gene_loss_type geneLossType)
// {
//   node->gene_loss_date_ = geneLossDate;
//   node->gene_loss_type_ = geneLossType;
//   node->protein_pointer_ = NULL;
//   for (int32_t i = 0; i < node->nb_promoters_; i++) {node->rna_pointers_[i] = NULL;}
//   if (geneLossDate > end_gener_) end_gener_ = geneLossDate;
//   nb_active_leaves_ --;
// }



GeneTreeNode *GeneTree::search_in_leaves(const Protein * protein)
{
  return root_->search_in_subtree_leaves(protein);
}


void GeneTree::print_to_screen()
{
  root_->print_subtree_to_screen();
  printf("\n");
}


void GeneTree::write_to_files(const char * topologyFileName, const char * nodeAttributesFileName, int32_t end_gener)
{
  FILE  * topology_file = fopen(topologyFileName, "w");
  if (topology_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file %s.\n", topologyFileName);
      exit(EXIT_FAILURE);
    }

  FILE  * node_attributes_file = fopen(nodeAttributesFileName, "w");
  if (node_attributes_file == NULL)
    {
      fprintf(stderr, "Error: cannot create file %s.\n", topologyFileName);
      exit(EXIT_FAILURE);
    }


  root_->write_subtree_to_files(topology_file, node_attributes_file, end_gener);
  fprintf(topology_file, ";");

  fclose(topology_file);
  fclose(node_attributes_file);
}


// f must already be open
void GeneTree::write_nodes_in_tabular_file(int32_t treeID, FILE * f)
{
  root_->write_subtree_nodes_in_tabular_file(treeID, f);
}


// =================================================================
//                           Protected Methods
// =================================================================

} // namespace aevol
