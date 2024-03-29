// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hydrate/RemoveFarAwayWatersMover.hh
/// @brief The RemoveWaters mover protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai, Lucas S. P. Rudden

#ifndef INCLUDED_protocols_hydrate_RemoveFarAwayWatersMover_HH
#define INCLUDED_protocols_hydrate_RemoveFarAwayWatersMover_HH

// Protocol headers
#include <protocols/hydrate/RemoveFarAwayWatersMover.fwd.hh>
#include <protocols/hydrate/RemoveWatersMover.hh>
#include <protocols/moves/Mover.hh>

// Core and Basic headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <iostream>
#include <string>

namespace protocols {
namespace hydrate {

/// @brief
class RemoveFarAwayWatersMover: public protocols::moves::Mover {
public:
	// Default constructor
	RemoveFarAwayWatersMover();

	// Copy constructor
	RemoveFarAwayWatersMover( RemoveFarAwayWatersMover const & hyd );

	// Assignment operator
	RemoveFarAwayWatersMover & operator=( RemoveFarAwayWatersMover const & hyd );

	// Destructor
	virtual ~RemoveFarAwayWatersMover();

	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;
	
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &
		) override;

	static void register_options();

	/// @brief  Apply the corresponding move to <input_pose>.
	virtual void apply(core::pose::Pose & input_pose);
	
	static std::string mover_name();

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_remove_far_away_water_mover( utility::tag::XMLSchemaDefinition & xsd );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
private:
	// Initialize data members.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(RemoveFarAwayWatersMover hyd_to, RemoveFarAwayWatersMover hyd_from);

};  // class RemoveFarAwayWatersMover

}  // namespace hydrate
}  // namespace protocols

#endif  // INCLUDED_protocols_hydrate_RemoveFarAwayWatersMover_HH
