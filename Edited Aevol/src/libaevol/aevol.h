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


#ifndef AEVOL_AEVOL_H_
#define AEVOL_AEVOL_H_

#include "ae_enums.h"
#include "ae_string.h"
#include "Utils.h"
#include "Alignment.h"
#include "Codon.h"
#include "Dna.h"
#include "DnaReplicationReport.h"
#include "Dump.h"
#include "ExpManager.h"
#include "ExpSetup.h"
#include "ParameterLine.h"
#include "Fuzzy.h"
#include "HybridFuzzy.h"
#include "Gaussian.h"
#include "GeneticUnit.h"
#include "GridCell.h"
#include "Habitat.h"
#include "IndividualFactory.h"
#include "Individual.h"
#include "Metrics.h"
#include "JumpingMT.h"
#include "JumpPoly.h"
#include "Logging.h"
#include "macros.h"
#include "Mutation.h"
#include "PointMutation.h"
#include "SmallInsertion.h"
#include "SmallDeletion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"
#include "MutationParams.h"
#include "NonCodingMetrics.h"
#include "OutputManager.h"
#include "ParamLoader.h"
#include "Phenotype.h"
#include "PhenotypicSegment.h"
#include "PhenotypicTarget.h"
#include "PhenotypicTargetHandler.h"
#include "Point.h"
#include "Protein.h"
#include "ReplicationReport.h"
#include "Rna.h"
#include "Selection.h"
#include "StatRecord.h"
#include "Stats.h"
#include "AeTime.h"
#include "Tree.h"
#include "VisAVis.h"
#include "World.h"

#include "LightTree.h"

#ifdef __X11
  #include "ExpManager_X11.h"
  #include "Individual_X11.h"
  #include "X11Window.h"
#endif

#endif // AEVOL_AEVOL_H_
