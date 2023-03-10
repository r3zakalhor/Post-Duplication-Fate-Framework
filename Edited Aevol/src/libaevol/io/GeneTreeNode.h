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

/** \class GeneTreeNode
 *  \brief Currently used only by post-treatments, on a specific lineage, to monitor the fate of paralogs.
 *         Each node corresponds to a coding RNA. When it is duplicated, two new nodes are added in the tree,
 *         as children of the ancestral version. The left child corresponds to the original DNA segment, while
 *         the right child corresponds to the copy that was reinserted elsewhere, possibly in another genetic
 *         unit.
 */


 #ifndef AEVOL_GENE_TREE_NODE_H_
#define  AEVOL_GENE_TREE_NODE_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>




// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Rna.h"
#include "GeneMutation.h"
#include "GeneTree.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================







class GeneTreeNode
{
  friend class GeneTree;

  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    GeneTreeNode(int32_t nodeCreationDate, Protein * protein);


    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~GeneTreeNode();

    // =================================================================
    //                            Public Methods
    // =================================================================
    GeneTreeNode * search_in_subtree_leaves(const Protein * protein);
    void print_subtree_to_screen(); // for debug purposes
    void write_subtree_to_files(FILE * topologyFile, FILE * nodeAttributesFile, int32_t end_gener);
    void write_subtree_nodes_in_tabular_file(int32_t treeID, FILE *f); // f must already be open
    void update_pointers_in_subtree_leaves(GeneticUnit * unit);
    void anticipate_mutation_effect_on_genes_in_subtree_leaves(const Mutation * mut, int32_t lengthOfGeneticUnit);
    void register_actual_mutation_effect_on_genes_in_subtree_leaves(GeneTree * tree, const Mutation * mut, GeneticUnit* unit, int32_t gener, double impact_on_metabolic_error);

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    GeneTreeNode()
      {
        printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      };
    GeneTreeNode(const GeneTreeNode &model)
      {
        printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      };

    // =================================================================
    //                          Protected Attributes
    // =================================================================
      static int32_t     nextID_;

    int32_t             ID_;
    int32_t             dna_creation_date_;   // generation when the DNA of this gene copy was created. For left nodes, equals the creation date of the parent node.
    int32_t             node_creation_date_;  // generation when this node was created. For right nodes, equals the dna creation date.
    int32_t             gene_loss_date_;      // generation when the gene was lost: became a pseudogene, or was deleted, or was broken by a rearrangement
    ae_gene_loss_type   gene_loss_type_;      // NOT_LOST_YET or LOST_BY_LOCAL_MUTATION or DELETED or BROKEN_BY_REAR
    Strand strand_;
    int32_t             shine_dal_position_;
    size_t             nb_promoters_;
    int32_t *           promoter_positions_;
    Protein *        protein_pointer_;     // for a leaf (current state of a gene), points to the potein object
                                              // for an internal node (ancestral state of a gene), points to NULL
                                              Rna **           rna_pointers_;         // for a leaf (current state of a gene), points to the RNA object
                                              // for an internal node (ancestral state of a gene), points to NULL
    std::list<GeneMutation *> mutation_list;       // list of ae_gene_mutations since the creation date of the node, i.e. since the last duplication


    GeneTreeNode * left_child_;    // NULL until the gene is duplicated, then points to the copy lying on the original DNA segment
    GeneTreeNode * right_child_;   // NULL until the gene is duplicated, then points to the copy lying on the duplicated DNA segment
                                        // (which was reinserted elsewhere in the genome, possibly on another genetic unit)
                                        GeneTreeNode * parent_node_;   // points to the node that corresponds to the state of the gene before the last duplication (NULL for the root of the gene tree)

    bool cds_possibly_modified_;
    bool cds_completely_deleted_;
    bool proms_possibly_modified_;
    bool gene_possibly_duplicated_;
    int32_t putative_position_for_the_duplicate_;


};
} // namespace aevol

#endif // AEVOL_GENE_TREE_NODE_H_
