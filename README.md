# EventGeneratorInterface
Start for an interface between neutrino generators: GENIE, NEUT, GiBUU, NuWro, ...

This essentially is based on neutgeom from NEUT, with some first work to include GiBUU.
Including GENIE should be fairly straightforward while for NuWro a separate function that 
calculates xsections is needed.

Ideally: pre-calculate xsec splines with different model variations for each generator and
interface the spline interpolation with TXSec.

Currently: event_rate.cc: Calculated total cross section and probabilities for each neutrino
	   		  flux vector in each part of the input geometry.
			  INPUTS: a TNeutrinoFlux
			  	  a TDetectorGeometry
				  a TXSec: generator specific. Given neutrino type, energy, A and Z, channel.
				    	   Interface for NEUT ready, first version for GiBUU in beta testing.
					   => Should all be spline lookups for speedup.
					   -> used by event_rate.cc through TNuTrajectory

			  OUTPUT : Rootfile with precalculated probabilities.

	  generate event (genev.cc): Use precalculate probabilities to choose interaction vertex and generate 
	  	   	 	     kinematics and final state particles.
				     INPUTS: TNeutrinoFlux
				     	     TDetectorGeometry
					     Rootfile from event_rate.cc
					     The type of output, eg. RooTracker, something flatter, a Les Houches, ... TMCOutput
				     uses TNuGenerator which interfaces with all generators to generate kinematics and
				     prepare output initial/final state particles for output format.


To Compile: same as standard neutgeom.