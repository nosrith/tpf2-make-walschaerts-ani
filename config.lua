return {
	-- basePoint: 
	--   14 reference points for the animation calculation.
	--   Specify the coordinates for each point at the neutral state.
	--
	--   You don't have to set all of these if you don't need animations 
	--   for all objects. The following table shows which reference points
	--   are needed for each object.
	--   
	--                     p1 p2 p3 p4 p5 p6 p7 p8 p9 10 11 12 13 14
	--   piston rod         R  R  R  R  -  -  -  -  -  -  -  -  -  -
	--   connection rod     R  R  R  R  -  -  -  -  -  -  -  -  -  -
	--   coupling rod       -  -  R  R  -  -  -  -  -  -  -  -  -  -
	--   return crank    -  -  -  R  R  -  -  -  -  -  -  -  -  -
	--   eccentric rod      -  -  -  R  R  R  R  -  -  -  -  -  -  -
	--   expansion link     -  -  -  R  R  R  R  -  -  -  -  -  -  -
	--   radius bar         R  R  R  R  R  R  R  R  R  R  R  R  -  -
	--   union link         R  R  R  R  R  R  R  R  R  R  R  R  -  -
	--   combination lever  R  R  R  R  R  R  R  R  R  R  R  R  -  -
	--   valve stem         R  R  R  R  R  R  R  R  R  R  R  R  -  -
	--   lifting link       R  R  R  R  R  R  R  R  R  R  R  R  R  R
	--
	--   Also don't have to set y values because they are not used at all.
	--
	basePoint = {
		-- p1: Center position of cylinder.
		p1 = { 3.485, 0.0, 0.8 },

		-- p2: Joint position of piston and connection rod.
		p2 = { 2.255, 0.0, 0.8 },

		-- p3: Joint position of connection and coupling rod (= crunk pin).
		p3 = { 0.0, 0.0, 0.495 },

		-- p4: Center position of main axle.
		p4 = { 0.0, 0.0, 0.8 },

		-- p5: Joint position of eccentric crank and eccentric rod.
		p5 = { 0.165, 0.0, 0.7808 },

		-- p6: Joint position of eccentric rod and expansion link.
		p6 = { 1.4836, 0.0, 0.9544 },

		-- p7: Axle position of expansion link.
		p7 = { 1.25, 0.0, 1.33 },

		-- p8: Joint position of expansion link and radius bar.
		-- p8 = { 1.3061, 0.0, 1.2197 },
		p8 = { 1.2635, 0.0, 1.2202 },

		-- p9: Joint position of crosshead arm and union link.
		p9 = { 2.255, 0.0, 0.5 },

		-- p10: Joint position of radius bar and combination lever.
		p10 = { 2.67, 0.0, 1.32 },

		-- p11: Joint position of union link and combination lever.
		p11 = { 2.62, 0.0, 0.5 },

		-- p12: Joint position of combination lever and valve stem.
		p12 = { 2.62, 0.0, 1.235 },

		-- p13: Joint position of radius bar and lifting link.
		p13 = { 1.5387, 0.0, 1.2535 },

		-- p14: Axle position of lifting link.
		p14 = { 1.49, 0.0, 1.6 },
	},

	-- origin: Origin coordinates for each node in global coordinate system.
	origin = {
		pistonRod = { 2.255, 0.0, 0.8 },
		connectionRod = { 0.0, 0.0, 0.495 },
		couplingRod = { 0.0, 0.0, 0.495 },
		returnCrank = { 0.0, 0.0, 0.495 },
		eccentricRod = { 0.0, 0.0, 0.825 },
		expansionLink = { 1.25, 0.0, 1.33 },
		radiusBar = { 0.0, 0.0, 0.0 },
		unionLink = { 2.255, 0.0, 0.8 },
		combinationLever = { 2.62, 0.0, 1.32 },
		valveStem = { 0.0, 0.0, 0.0 },
		liftingLink = { 1.49, 0.0, 1.6 },
	},

	-- phaseOffset: 
	--   Wheel phase offset at t0. 0.0 means no offset and 1.0 and -1.0 means completely rounded.
	--   For example set 0.0 for right side ani, and set -0.25 for left to delay 90 degrees for cross balancing.
	phaseOffset = 0.0,

	-- nstep: The number of steps the wheels take around.
	nstep = 16,

	-- outputDirPath: If empty files will be generated at current directry. 
	outputDirPath = "",

	-- outputName: If empty files will not be generated.
	outputName = {
		pistonRod = "piston_rod",
		connectionRod = "connection_rod",
		couplingRod = "coupling_rod",
		returnCrank = "return_crank",
		eccentricRod = "eccentric_rod",
		expansionLink = "expansion_link",
		radiusBar = "radius_bar",
		unionLink = "union_link",
		combinationLever = "combination_lever",
		valveStem = "valve_stem",
		liftingLink = "lifting_link",
	},

	-- outputSuffix: Useful to distinguish between right and left.
	outputSuffix = ".ani",

	-- verbose: Whether show details on console.
	verbose = true,
}
