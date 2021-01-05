local lbfgs = require "lbfgs"
local vec2 = require "vec2"

local config
local ctx = {}

local zero2 = vec2.new(0.0, 0.0)

local function sq(x)
	return x * x
end

local function checkVars(...)
	for _, v in ipairs({...}) do
		if not ctx[v] then return false end
	end
	return true
end

local function log(s, ...)
	if config.verbose then
		print(string.format(s, ...))
	end
end

local function calcP3()
	if not checkVars("bp3", "bp4") then 
		ctx.p3 = nil
		return
	end
	ctx.p3 = (ctx.bp3 - ctx.bp4):rot(-(ctx.tfrac + config.phaseOffset) * 2 * math.pi) + ctx.bp4
	log("  P3: %s (BP3: %s, d: %s)", ctx.p3, ctx.bp3, ctx.p3 - ctx.bp3)
end

local function calcP2()
	if not checkVars("bp1", "bp2", "bp3", "p3") then 
		ctx.p2 = nil
		return 
	end
	local bp21 = ctx.bp1 - ctx.bp2
	local blen21 = bp21:len()
	local sqblen32 = (ctx.bp2 - ctx.bp3):len2()
	local d = bp21:dot(ctx.bp2 - ctx.p3)
	local fortho = -d / sq(blen21)
	local fpara = math.sqrt(sqblen32 - (ctx.bp2 - ctx.p3):len2() + sq(d / blen21)) / blen21
	ctx.p2 = ctx.bp2 + bp21 * (fortho + fpara)
	log("  P2: %s (BP2: %s, d: %s)", ctx.p2, ctx.bp2, ctx.p2 - ctx.bp2)
end

local function calcP1()
	if not checkVars("bp1", "bp2", "p2") then
		ctx.p1 = nil
		return
	end
	ctx.p1 = ctx.p2 + (ctx.bp1 - ctx.bp2)
	log("  P1: %s (BP1: %s, d: %s)", ctx.p1, ctx.bp1, ctx.p1 - ctx.bp1)
end

local function calcP5()
	if not checkVars("bp4", "bp5") then
		ctx.p5 = nil
		return
	end
	ctx.p5 = (ctx.bp5 - ctx.bp4):rot(-(ctx.tfrac + config.phaseOffset) * 2 * math.pi) + ctx.bp4
	log("  P5: %s (BP5: %s, d: %s)", ctx.p5, ctx.bp5, ctx.p5 - ctx.bp5)
end

local function calcP6()
	if not checkVars("bp5", "bp6", "bp7", "p5", "p6") then
		ctx.p6 = nil
		return
	end
	local sqblen56 = (ctx.bp6 - ctx.bp5):len2()
	local sqblen76 = (ctx.bp6 - ctx.bp7):len2()
	local eval = function(x, g, n, step)
		local p6 = vec2.new(x[1], x[2])
		local sqlen56 = (p6 - ctx.p5):len2()
		local sqlen76 = (p6 - ctx.bp7):len2()
		local c1 = sqlen56 - sqblen56
		local c2 = sqlen76 - sqblen76
		local dp6 = ((p6 - ctx.p5) * c1 + (p6 - ctx.bp7) * c2) * 4
		g[1], g[2] = dp6[1], dp6[2]
		return sq(c1) + sq(c2)
	end
	local ret, x = lbfgs.lbfgs({ ctx.p6[1], ctx.p6[2] }, eval)
	if ret < 0 then
		print("**ERROR**")
		print("Calculation failed at P6 at step " .. ctx.t .. ".")
		print("This often occures when the base points are in an inoperable position.")
		print("Internal message: ")
		print("  " .. ret .. ": " .. lbfgs.strerror(ret))
		ctx.error = true
		ctx.p6 = nil
		return
	end
	ctx.p6 = vec2.new(x[1], x[2])
	log("  P6: %s (BP6: %s, d: %s)", ctx.p6, ctx.bp6, ctx.p6 - ctx.bp6)
end

local function calcP8()
	if not checkVars("bp7", "p6") then
		ctx.p8 = nil
		return
	end
	local bang6 = (ctx.bp6 - ctx.bp7):angle()
	local ang6 = (ctx.p6 - ctx.bp7):angle()
	ctx.p8 = (ctx.bp8 - ctx.bp7):rot(ang6 - bang6) + ctx.bp7
	log("  P8: %s (BP8: %s, d: %s)", ctx.p8, ctx.bp8, ctx.p8 - ctx.bp8)
end

local function calcP9()
	if not checkVars("bp2", "bp9", "p2") then
		ctx.p9 = nil
		return
	end
	ctx.p9 = ctx.p2 + (ctx.bp9 - ctx.bp2)
	log("  P9: %s (BP9: %s, d: %s)", ctx.p9, ctx.bp9, ctx.p9 - ctx.bp9)
end

local function calcP10P11P12()
	if not checkVars("bp1", "bp2", "bp8", "bp9", "bp10", "bp11", "bp12", "p8", "p9", "p10", "p11") then
		ctx.p10, ctx.p11, ctx.p12 = nil, nil, nil
		return
	end
	local sqblen810 = (ctx.bp10 - ctx.bp8):len2()
	local sqblen911 = (ctx.bp11 - ctx.bp9):len2()
	local sqblen1011 = (ctx.bp11 - ctx.bp10):len2()
	local bp21 = ctx.bp1 - ctx.bp2
	local eval = function(x, g, n, step)
		local p10 = vec2.new(x[1], x[2])
		local p11 = vec2.new(x[3], x[4])
		local f12to21 = x[5]
		local sqlen810 = (p10 - ctx.p8):len2()
		local sqlen911 = (p11 - ctx.p9):len2()
		local sqlen1011 = (p11 - p10):len2()
		local p12by1011 = p10 + (p11 - p10) * ctx.f12to1011
		local p12by21 = ctx.bp12 + bp21 * f12to21
		local c1 = sqlen810 - sqblen810
		local c2 = sqlen911 - sqblen911
		local c3 = sqlen1011 - sqblen1011
		local c4 = p12by1011 - p12by21
		local dp10 = ((p10 - ctx.p8) * c1 - (p11 - p10) * c3) * 4 + c4 * (1 - ctx.f12to1011) * 2
		local dp11 = ((p11 - ctx.p9) * c2 + (p11 - p10) * c3) * 4 + c4 * ctx.f12to1011 * 2
		local df12to21 = -2 * c4:dot(bp21)
		g[1], g[2], g[3], g[4], g[5] = dp10[1], dp10[2], dp11[1], dp11[2], df12to21
		return sq(c1) + sq(c2) + sq(c3) + c4:len2()
	end
	local ret, x = lbfgs.lbfgs({ ctx.p10[1], ctx.p10[2], ctx.p11[1], ctx.p11[2], (ctx.p12 - ctx.bp12):len() / bp21:len() }, eval)
	if ret < 0 then
		print("**ERROR**")
		print("Calculation failed at P10, P11 and P12 at step " .. ctx.t .. ".")
		print("This often occures when the base points are in an inoperable position.")
		print("Internal message: ")
		print("  " .. ret .. ": " .. lbfgs.strerror(ret))
		ctx.error = true
		ctx.p10, ctx.p11, ctx.p12 = nil, nil, nil
		return
	end
	ctx.p10 = vec2.new(x[1], x[2])
	ctx.p11 = vec2.new(x[3], x[4])
	ctx.p12 = ctx.p10 + (ctx.p11 - ctx.p10) * ctx.f12to1011
	log("  P10: %s (BP10: %s, d: %s)", ctx.p10, ctx.bp10, ctx.p10 - ctx.bp10)
	log("  P11: %s (BP11: %s, d: %s)", ctx.p11, ctx.bp11, ctx.p11 - ctx.bp11)
	log("  P12: %s (BP12: %s, d: %s)", ctx.p12, ctx.bp12, ctx.p12 - ctx.bp12)
end

local function calcP13()
	if not checkVars("p8", "p10") then
		ctx.p13 = nil
		return
	end
	ctx.p13 = ctx.p8 + (ctx.p10 - ctx.p8) * ctx.f13to810
	log("  P13: %s (BP13: %s, d: %s)", ctx.p13, ctx.bp13, ctx.p13 - ctx.bp13)
end

local function makeTransf(rotY, cnt, transl, name)
	log("  %s: {rotY = %f, cnt = %s, transl = %s}", name, rotY, cnt, transl)

	local sy = math.sin(rotY)
	local cy = math.cos(rotY)

	local tx = cnt[1] + transl[1]
	local tz = cnt[2] + transl[2]

	return {
		cy, 0.0, sy, 0.0,
		0.0, 1.0, 0.0, 0.0,
		-sy, 0.0, cy, 0.0,
		-cy * cnt[1] + sy * cnt[2] + tx, 0.0, - sy * cnt[1] - cy * cnt[2] + tz, 1.0
	}
end

local function calcPistonRodTransf()
	if not (ctx.transf.pistonRod and checkVars("p2")) then
		ctx.transf.pistonRod = nil
		return
	end
	ctx.transf.pistonRod[ctx.t + 1] = makeTransf(0.0, zero2, ctx.p2 - ctx.bp2, "piston rod")
end

local function calcConnectingRodTransf()
	if not (ctx.transf.connectingRod and checkVars("p2", "p3")) then
		ctx.transf.connectingRod = nil
		return
	end
	ctx.transf.connectingRod[ctx.t + 1] = makeTransf(
		(ctx.p2 - ctx.p3):angle() - (ctx.bp2 - ctx.bp3):angle(),
		ctx.bp3 - vec2.fromVec3(config.origin.connectingRod),
		ctx.p3 - ctx.bp3,
		"connecting rod"
	)
end

local function calcCouplingRodTransf()
	if not (ctx.transf.couplingRod and checkVars("p3")) then
		ctx.transf.couplingRod = nil
		return
	end
	ctx.transf.couplingRod[ctx.t + 1] = makeTransf(0.0, zero2, ctx.p3 - ctx.bp3, "coupling rod")
end

local function calcReturnCrankTransf()
	if not (ctx.transf.returnCrank and checkVars("p3", "p5")) then
		ctx.transf.returnCrank = nil
		return
	end
	ctx.transf.returnCrank[ctx.t + 1] = makeTransf(
		(ctx.p5 - ctx.p3):angle() - (ctx.bp5 - ctx.bp3):angle(),
		ctx.bp3 - vec2.fromVec3(config.origin.returnCrank),
		ctx.p3 - ctx.bp3,
		"return crank"
	)
end

local function calcEccentricRodTransf()
	if not (ctx.transf.eccentricRod and checkVars("p5", "p6")) then
		ctx.transf.eccentricRod = nil
		return
	end
	ctx.transf.eccentricRod[ctx.t + 1] = makeTransf(
		(ctx.p6 - ctx.p5):angle() - (ctx.bp6 - ctx.bp5):angle(),
		ctx.bp5 - vec2.fromVec3(config.origin.eccentricRod),
		ctx.p5 - ctx.bp5,
		"eccentric rod"
	)
end

local function calcExpansionLinkTransf()
	if not (ctx.transf.expansionLink and checkVars("p6")) then
		ctx.transf.expansionLink = nil
		return
	end
	ctx.transf.expansionLink[ctx.t + 1] = makeTransf(
		(ctx.p6 - ctx.bp7):angle() - (ctx.bp6 - ctx.bp7):angle(),
		ctx.bp7 - vec2.fromVec3(config.origin.expansionLink),
		zero2,
		"expansion link"
	)
end

local function calcRadiusBarTransf()
	if not (ctx.transf.radiusBar and checkVars("p8", "p10")) then
		ctx.transf.radiusBar = nil
		return
	end
	ctx.transf.radiusBar[ctx.t + 1] = makeTransf(
		(ctx.p10 - ctx.p8):angle() - (ctx.bp10 - ctx.bp8):angle(),
		ctx.bp8 - vec2.fromVec3(config.origin.radiusBar),
		ctx.p8 - ctx.bp8,
		"radius bar"
	)
end

local function calcUnionLinkTransf()
	if not (ctx.transf.unionLink and checkVars("p9", "p11")) then
		ctx.transf.unionLink = nil
		return
	end
	ctx.transf.unionLink[ctx.t + 1] = makeTransf(
		(ctx.p11 - ctx.p9):angle() - (ctx.bp11 - ctx.bp9):angle(),
		ctx.bp9 - vec2.fromVec3(config.origin.unionLink),
		ctx.p9 - ctx.bp9,
		"union link"
	)
end

local function calcCombinationLeverTransf()
	if not (ctx.transf.combinationLever and checkVars("p10", "p11")) then
		ctx.transf.combinationLever = nil
		return
	end
	ctx.transf.combinationLever[ctx.t + 1] = makeTransf(
		(ctx.p11 - ctx.p10):angle() - (ctx.bp11 - ctx.bp10):angle(),
		ctx.bp10 - vec2.fromVec3(config.origin.combinationLever),
		ctx.p10 - ctx.bp10,
		"combination lever"
	)
end

local function calcValveStemTransf()
	if not (ctx.transf.valveStem and checkVars("p12")) then
		ctx.transf.valveStem = nil
		return
	end
	ctx.transf.valveStem[ctx.t + 1] = makeTransf(0.0, zero2, ctx.p12 - ctx.bp12, "valve stem")
end

local function calcLiftingLinkTransf()
	if not (ctx.transf.liftingLink and checkVars("p13")) then
		ctx.transf.liftingLink = nil
		return
	end
	ctx.transf.liftingLink[ctx.t + 1] = makeTransf(
		(ctx.p13 - ctx.bp14):angle() - (ctx.bp13 - ctx.bp14):angle(),
		ctx.bp14 - vec2.fromVec3(config.origin.liftingLink),
		zero2,
		"lifting link"
	)
end

local function step(t)
	log("[Step %d]", t)
	ctx.t = t
	ctx.tfrac = t / config.nstep
	ctx.times[t + 1] = ctx.tfrac * ctx.duration

	calcP3()
	calcP2()
	calcP1()
	calcP5()
	calcP6()
	calcP8()
	calcP9()
	calcP10P11P12()
	calcP13()

	calcPistonRodTransf()
	calcConnectingRodTransf()
	calcCouplingRodTransf()
	calcReturnCrankTransf()
	calcEccentricRodTransf()
	calcExpansionLinkTransf()
	calcRadiusBarTransf()
	calcUnionLinkTransf()
	calcCombinationLeverTransf()
	calcValveStemTransf()
	calcLiftingLinkTransf()
end

local function setupContext()
	log("[Preprocessing]")

	-- Make vec2 instance for all base points
	for k, v in pairs(config.basePoint) do
		ctx["b" .. k] = v and vec2.fromVec3(v) or nil
	end

	-- Align p12 onto p10-p11 line
	if ctx.bp1 and ctx.bp2 and ctx.bp10 and ctx.bp11 and ctx.bp12 then
		local bp120 = ctx.bp12
		local bp1011 = ctx.bp11 - ctx.bp10
		local bp1012 = ctx.bp12 - ctx.bp10
		local bp21 = ctx.bp1 - ctx.bp2
		ctx.f12to1011 =
			(bp1012[1] * bp21[2] - bp1012[2] * bp21[1]) /
			(bp1011[1] * bp21[2] - bp1011[2] * bp21[1])
		ctx.bp12 = ctx.bp10 + bp1011 * ctx.f12to1011
		log("  Adjusted BP12: %s (BP12: %s, d: %s)", ctx.bp12, bp120, ctx.bp12 - bp120)
	end

	-- Align p13 onto p8-p10 line
	if ctx.bp8 and ctx.bp10 and ctx.bp13 and ctx.bp14 then
		local bp130 = ctx.bp13
		local bp810 = ctx.bp10 - ctx.bp8
		local bp814 = ctx.bp14 - ctx.bp8
		local bp1413 = ctx.bp13 - ctx.bp14
		ctx.f13to810 = 
			(bp814[1] * bp1413[2] - bp814[2] * bp1413[1]) /
			(bp810[1] * bp1413[2] - bp810[2] * bp1413[1])
		ctx.bp13 = ctx.bp8 + bp810 * ctx.f13to810
		log("  Adjusted BP13: %s (BP13: %s, d: %s)", ctx.bp13, bp130, ctx.bp13 - bp130)
	end

	-- set initial value for L-BFGS (and update by step)
	ctx.p6 = ctx.bp6
	ctx.p10 = ctx.bp10
	ctx.p11 = ctx.bp11
	ctx.p12 = ctx.bp12

	ctx.transf = {
		pistonRod = {},
		connectingRod = {},
		couplingRod = {},
		returnCrank = {},
		eccentricRod = {},
		expansionLink = {},
		radiusBar = {},
		unionLink = {},
		combinationLever = {},
		valveStem = {},
		liftingLink = {},
	}

	ctx.duration = 2 * math.pi * ctx.bp4[2] * 1000
	ctx.times = {}
end

local function makeAniFile(k)
	if not ctx.transf[k] or not config.outputName[k] then
		return
	end

	local prefix = config.outputDirPath and config.outputDirPath ~= "" and config.outputDirPath .. "/" or ""
	local filepath = prefix .. config.outputName[k] .. config.outputSuffix
	local f = assert(io.open(filepath, "w"))
	f:write("function data()\nreturn {\n    times = { ")
	f:write(table.concat(ctx.times, ", "))
	f:write(" },\n        transfs = {\n")
	for _, m in ipairs(ctx.transf[k]) do
		f:write("        { ")
		f:write(table.concat(m, ", "))
		f:write(" },\n")
	end
	f:write("    },\n}\nend\n")
	f:close()
end

local function main()
	config = dofile(arg[1] or "config.lua")

	setupContext()
	for t = 0, config.nstep do
		step(t)
	end

	for k in pairs(config.outputName) do
		makeAniFile(k)
	end

	return not ctx.error
end

local stat, msg = xpcall(main, function(err)
	print("**ERROR**")
	print(err)
	print(debug.traceback())
end)
if not stat or not msg then
	print("An error has occured. Check logs for detail.")
	io.read()
end
