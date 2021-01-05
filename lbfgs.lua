--
-- This is a Lua library of L-BFGS reimplemented from liblbfgs:
-- https://github.com/chokkan/liblbfgs
--
--
--      C library of Limited memory BFGS (L-BFGS).
--
-- Copyright (c) 1990, Jorge Nocedal
-- Copyright (c) 2007-2010 Naoaki Okazaki
-- All rights reserved.
--
-- Permission is hereby granted, free of charge, to any person obtaining a copy
-- of this software and associated documentation files (the "Software"), to deal
-- in the Software without restriction, including without limitation the rights
-- to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
-- copies of the Software, and to permit persons to whom the Software is
-- furnished to do so, subject to the following conditions:
--
-- The above copyright notice and this permission notice shall be included in
-- all copies or substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
-- IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
-- FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
-- AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
-- LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
-- OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
-- THE SOFTWARE.
--

local lib = {}

-- Return values of lbfgs().
-- Roughly speaking, a negative value indicates an error.
lib.SUCCESS = 0
lib.CONVERGENCE = 0
lib.STOP = 1
lib.ALREADY_MINIMIZED = 2
lib.ERR_UNKNOWNERROR = -1024
lib.ERR_LOGICERROR = -1023
lib.ERR_OUTOFMEMORY = -1022
lib.ERR_CANCELED = -1021
lib.ERR_INVALID_N = -1020
lib.ERR_INVALID_N_SSE = -1019
lib.ERR_INVALID_X_SSE = -1018
lib.ERR_INVALID_EPSILON = -1017
lib.ERR_INVALID_TESTPERIOD = -1016
lib.ERR_INVALID_DELTA = -1015
lib.ERR_INVALID_LINESEARCH = -1014
lib.ERR_INVALID_MINSTEP = -1013
lib.ERR_INVALID_MAXSTEP = -1012
lib.ERR_INVALID_FTOL = -1011
lib.ERR_INVALID_WOLFE = -1010
lib.ERR_INVALID_GTOL = -1009
lib.ERR_INVALID_XTOL = -1008
lib.ERR_INVALID_MAXLINESEARCH = -1007
lib.ERR_INVALID_ORTHANTWISE = -1006
lib.ERR_INVALID_ORTHANTWISE_START = -1005
lib.ERR_INVALID_ORTHANTWISE_END = -1004
lib.ERR_OUTOFINTERVAL = -1003
lib.ERR_INCORRECT_TMINMAX = -1002
lib.ERR_ROUNDING_ERROR = -1001
lib.ERR_MINIMUMSTEP = -1000
lib.ERR_MAXIMUMSTEP = -999
lib.ERR_MAXIMUMLINESEARCH = -998
lib.ERR_MAXIMUMITERATION = -997
lib.ERR_WIDTHTOOSMALL = -996
lib.ERR_INVALIDPARAMETERS = -995
lib.ERR_INCREASEGRADIENT = -994

-- Line search algorithms.
lib.LINESEARCH_DEFAULT = 0
lib.LINESEARCH_MORETHUENTE = 0
lib.LINESEARCH_BACKTRACKING_ARMIJO = 1
lib.LINESEARCH_BACKTRACKING = 2
lib.LINESEARCH_BACKTRACKING_WOLFE = 2
lib.LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3

local defparam = {
	m = 6,
	epsilon = 1e-5,
	past = 0,
	delta = 1e-5,
	max_iterations = 0,
	linesearch = 0,
	max_linesearch = 40,
	min_step = 1e-20,
	max_step = 1e20,
	ftol = 1e-4,
	wolfe = 0.9,
	gtol = 0.9,
	xtol = 1e-16,
	orthantwise_c = 0.0,
	orthantwise_start = 1,
	orthantwise_end = -1,
}

local function fsigndiff(x, y)
	return x * (y / math.abs(y)) < 0.
end

local function vecset(x, c, n)
	for i = 1, n do 
		x[i] = c 
	end
end

local function veccpy(y, x, n)
	for i = 1, n do
		y[i] = x[i]
	end
end

local function vecncpy(y, x, n)
	for i = 1, n do
		y[i] = -x[i]
	end
end

local function vecadd(y, x, c, n)
	for i = 1, n do
		y[i] = y[i] + c * x[i]
	end
end

local function vecdiff(z, x, y, n)
	for i = 1, n do
		z[i] = x[i] - y[i]
	end
end

local function vecscale(y, c, n)
	for i = 1, n do
		y[i] = y[i] * c
	end
end

local function vecmul(y, x, n)
	for i = 1, n do
		y[i] = y[i] * x[i]
	end
end

local function vecdot(x, y, n)
	local s = 0.
	for i = 1, n do
		s = s + x[i] * y[i]
	end
	return s
end

local function vec2norm(x, n)
	return math.sqrt(vecdot(x, x, n))
end

local function vec2norminv(x, n)
	return 1.0 / vec2norm(x, n)
end

local function owlqn_x1norm(x, start, n)
	local norm = 0.
	for i = start, n do
		norm = norm + math.abs(x[i])
	end
	return norm
end

local function owlqn_pseudo_gradient(pg, x, g, n, c, start, end_)
	-- Compute the negative of gradients.
	for i = 1, start - 1 do
		pg[i] = g[i]
	end

	-- Compute the psuedo-gradients.
	for i = start, end_ do
		if x[i] < 0. then
			-- Differentiable.
			pg[i] = g[i] - c
		elseif 0. < x[i] then
			-- Differentiable.
			pg[i] = g[i] + c
		else
			if g[i] < -c then
				-- Take the right partial derivative.
				pg[i] = g[i] + c
			elseif c < g[i] then
				-- Take the left partial derivative.
				pg[i] = g[i] - c
			else
				pg[i] = 0.
			end
		end
	end

	for i = end_ + 1, n do
		pg[i] = g[i]
	end
end

local function owlqn_project(d, sign, start, end_)
	for i = start, end_ do
        if d[i] * sign[i] <= 0 then
            d[i] = 0
		end
    end
end

local function line_search_backtracking(n, x, f, g, s, stp, xp, gp, wp, evaluate, param)
	local count = 0
	local dec = 0.5
	local inc = 2.1

	-- Check the input parameters for errors.
	if stp <= 0. then
		return lib.ERR_INVALIDPARAMETERS, f, stp
	end

    -- Compute the initial gradient in the search direction.
    local dginit = vecdot(g, s, n)

    -- Make sure that s points to a descent direction.
    if 0 < dginit then
        return lib.ERR_INCREASEGRADIENT, f, stp
	end

    -- The initial value of the objective function.
    local finit = f
    local dgtest = param.ftol * dginit
    while true do
        veccpy(x, xp, n)
        vecadd(x, s, stp, n)

        -- Evaluate the function and gradient values.
        f = evaluate(x, g, n, stp)

        count = count + 1

		local width
        if f > finit + stp * dgtest then
            width = dec
		else
            -- The sufficient decrease condition (Armijo condition).
            if param.linesearch == lib.LINESEARCH_BACKTRACKING_ARMIJO then
                -- Exit with the Armijo condition.
                return count, f, stp
			end

	        -- Check the Wolfe condition.
	        local dg = vecdot(g, s, n)
	        if dg < param.wolfe * dginit then
    		    width = inc
			else
		        if param.linesearch == lib.LINESEARCH_BACKTRACKING_WOLFE then
		            -- Exit with the regular Wolfe condition.
		            return count, f, stp
				end

		        -- Check the strong Wolfe condition.
		        if dg > -param.wolfe * dginit then
		            width = dec
				else
		            -- Exit with the strong Wolfe condition.
		            return count, f, stp
		        end
            end
        end

        if stp < param.min_step then
			-- The step is the minimum value.
            return lib.ERR_MINIMUMSTEP
		end
        if stp > param.max_step then
            -- The step is the maximum value.
            return lib.ERR_MAXIMUMSTEP
		end
        if param.max_linesearch <= count then
            -- Maximum number of iteration.
            return lib.ERR_MAXIMUMLINESEARCH
		end

		stp = stp * width
    end
end

local function line_search_backtracking_owlqn(n, x, f, g, s, stp, xp, gp, wp, evaluate, param)
	local count = 0
	local width = 0.5
	local finit = f

	-- Check the input parameters for errors.
	if stp <= 0. then
		return lib.ERR_INVALIDPARAMETERS, f, stp
	end

	-- Choose the orthant for the new point.
	for i = 1, n do
        wp[i] = xp[i] == 0. and -gp[i] or xp[i]
    end

	while true do
        -- Update the current point.
        veccpy(x, xp, n)
        vecadd(x, s, stp, n)

        -- The current point is projected onto the orthant.
        owlqn_project(x, wp, param.orthantwise_start, param.orthantwise_end)

        -- Evaluate the function and gradient values.
        f = evaluate(x, g, n, stp)

        -- Compute the L1 norm of the variables and add it to the object value.
        local norm = owlqn_x1norm(x, param.orthantwise_start, param.orthantwise_end)
        f = f + norm * param.orthantwise_c

        count = count + 1

        local dgtest = 0.
        for i = 1, n do
            dgtest = dgtest + (x[i] - xp[i]) * gp[i]
		end

        if f <= finit + param.ftol * dgtest then
            -- The sufficient decrease condition.
            return count, f, stp
		end

        if stp < param.min_step then
			-- The step is the minimum value.
            return lib.ERR_MINIMUMSTEP
		end
        if stp > param.max_step then
            -- The step is the maximum value.
            return lib.ERR_MAXIMUMSTEP
		end
        if param.max_linesearch <= count then
            -- Maximum number of iteration.
            return lib.ERR_MAXIMUMLINESEARCH
		end

        stp = stp * width
    end
end

---
-- Find a minimizer of an interpolated cubic function.
--  @param  u       The value of one point, u.
--  @param  fu      The value of f(u).
--  @param  du      The value of f'(u).
--  @param  v       The value of another point, v.
--  @param  fv      The value of f(v).
--  @param  du      The value of f'(v).
--  @return         The minimizer of the interpolated cubic.
local function cubic_minimizer(u, fu, du, v, fv, dv)
    local d = v - u
    local theta = (fu - fv) * 3 / d + du + dv
    local p = math.abs(theta)
    local q = math.abs(du)
    local r = math.abs(dv)
    local s = math.max(p, q, r)
    local a = theta / s
    local gamma = s * math.sqrt(a * a - (du / s) * (dv / s))
    if v < u then gamma = -gamma end
    p = gamma - du + theta
    q = gamma - du + gamma + dv
    r = p / q
	cm = u + r * d
	return cm
end

---
-- Find a minimizer of an interpolated cubic function.
--  @param  u       The value of one point, u.
--  @param  fu      The value of f(u).
--  @param  du      The value of f'(u).
--  @param  v       The value of another point, v.
--  @param  fv      The value of f(v).
--  @param  du      The value of f'(v).
--  @param  xmin    The maximum value.
--  @param  xmin    The minimum value.
--  @return         The minimizer of the interpolated cubic.
local function cubic_minimizer2(u, fu, du, v, fv, dv, xmin, xmax)
    local d = v - u
    local theta = (fu - fv) * 3 / d + du + dv
    local p = math.abs(theta)
    local q = math.abs(du)
    local r = math.abs(dv)
    local s = math.max(p, q, r)
    local a = theta / s
    local gamma = s * math.sqrt(math.max(0, a * a - (du / s) * (dv / s)))
    if u < v then gamma = -gamma end
    p = gamma - dv + theta
    q = gamma - dv + gamma + du
	r = p / q
	local cm
    if r < 0. and gamma ~= 0. then
        cm = v - r * d
    elseif a < 0 then
        cm = xmax
    else
        cm = xmin
	end
	return cm
end

---
-- Find a minimizer of an interpolated quadratic function.
--  @param  u       The value of one point, u.
--  @param  fu      The value of f(u).
--  @param  du      The value of f'(u).
--  @param  v       The value of another point, v.
--  @param  fv      The value of f(v).
--  @return         The minimizer of the interpolated quadratic.
local function quard_minimizer(u, fu, du, v, fv)
    local a = v - u
	local qm = u + du / ((fu - fv) / a + du) / 2 * a
	return qm
end

---
-- Find a minimizer of an interpolated quadratic function.
--  @param  u       The value of one point, u.
--  @param  du      The value of f'(u).
--  @param  v       The value of another point, v.
--  @param  dv      The value of f'(v).
--  @return         The minimizer of the interpolated quadratic.
local function quard_minimizer2(u, du, v, dv)
    local a = u - v
	local qm = v + dv / (dv - du) * a
	return qm
end

---
-- Update a safeguarded trial value and interval for line search.
--  The parameter x represents the step with the least function value.
--  The parameter t represents the current step. This function assumes
--  that the derivative at the point of x in the direction of the step.
--  If the bracket is set to true, the minimizer has been bracketed in
--  an interval of uncertainty with endpoints between x and y.
--  @param  x       The pointer to the value of one endpoint.
--  @param  fx      The pointer to the value of f(x).
--  @param  dx      The pointer to the value of f'(x).
--  @param  y       The pointer to the value of another endpoint.
--  @param  fy      The pointer to the value of f(y).
--  @param  dy      The pointer to the value of f'(y).
--  @param  t       The pointer to the value of the trial value, t.
--  @param  ft      The pointer to the value of f(t).
--  @param  dt      The pointer to the value of f'(t).
--  @param  tmin    The minimum value for the trial value, t.
--  @param  tmax    The maximum value for the trial value, t.
--  @param  brackt  The pointer to the predicate if the trial value is bracketed.
--  @return int     Status value. Zero indicates a normal termination.
--  
--  @see
--      Jorge J. More and David J. Thuente. Line search algorithm with
--      guaranteed sufficient decrease. ACM Transactions on Mathematical
--      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
local function update_trial_interval(x, fx, dx, y, fy, dy, t, ft, dt, tmin, tmax, brackt)
	local dsign = fsigndiff(dt, dx)

    -- Check the input parameters for errors.
    if brackt ~= 0 then
        if t <= math.min(x, y) or math.max(x, y) <= t then
            -- The trival value t is out of the interval.
            return lib.ERR_OUTOFINTERVAL, x, fx, dx, y, fy, dy, t, ft, dt, brackt
		end
		if 0. <= dx * (t - x) then
            -- The function must decrease from x.
            return lib.ERR_INCREASEGRADIENT, x, fx, dx, y, fy, dy, t, ft, dt, brackt
		end
        if tmax < tmin then
            -- Incorrect tmin and tmax specified.
            return lib.ERR_INCORRECT_TMINMAX, x, fx, dx, y, fy, dy, t, ft, dt, brackt
		end
    end

	-- Trial value selection.
	local bound
	local mc  -- minimizer of an interpolated cubic.
	local mq  -- minimizer of an interpolated quadratic.
	local newt  -- new trial value.
    if fx < ft then
		-- Case 1: a higher function value.
		-- The minimum is brackt. If the cubic minimizer is closer
		-- to x than the quadratic one, the cubic one is taken, else
		-- the average of the minimizers is taken.
        brackt = 1
        bound = 1
        mc = cubic_minimizer(x, fx, dx, t, ft, dt)
        mq = quard_minimizer(x, fx, dx, t, ft)
        if math.abs(mc - x) < math.abs(mq - x) then
			newt = mc
		else
            newt = mc + 0.5 * (mq - mc)
        end
    elseif dsign ~= 0 then
		-- Case 2: a lower function value and derivatives of
		-- opposite sign. The minimum is brackt. If the cubic
		-- minimizer is closer to x than the quadratic (secant) one,
		-- the cubic one is taken, else the quadratic one is taken.
        brackt = 1
		bound = 0
		mc = cubic_minimizer(x, fx, dx, t, ft, dt)
        mq = quard_minimizer2(x, dx, t, dt)
        if math.abs(mc - t) > math.abs(mq - t) then
            newt = mc
		else
            newt = mq
        end
    elseif math.abs(dt) < math.abs(dx) then
		-- Case 3: a lower function value, derivatives of the
		-- same sign, and the magnitude of the derivative decreases.
		-- The cubic minimizer is only used if the cubic tends to
		-- infinity in the direction of the minimizer or if the minimum
		-- of the cubic is beyond t. Otherwise the cubic minimizer is
		-- defined to be either tmin or tmax. The quadratic (secant)
		-- minimizer is also computed and if the minimum is brackt
		-- then the the minimizer closest to x is taken, else the one
		-- farthest away is taken.
        bound = 1
        mc = cubic_minimizer2(x, fx, dx, t, ft, dt, tmin, tmax)
        mq = quard_minimizer2(x, dx, t, dt)
        if brackt ~= 0 then
            if math.abs(t - mc) < math.abs(t - mq) then
                newt = mc
			else
                newt = mq;
			end
		else
            if math.abs(t - mc) > math.abs(t - mq) then
                newt = mc
			else
                newt = mq;
            end
		end
	else
		-- Case 4: a lower function value, derivatives of the
		-- same sign, and the magnitude of the derivative does
		-- not decrease. If the minimum is not brackt, the step
		-- is either tmin or tmax, else the cubic minimizer is taken.
        bound = 0
        if brackt ~= 0 then
            newt = cubic_minimizer(t, ft, dt, y, fy, dy)
        elseif x < t then
            newt = tmax
		else
            newt = tmin
        end
    end

	-- Update the interval of uncertainty. This update does not
	-- depend on the new step or the case analysis above.

	-- - Case a: if f(x) < f(t),
	--	 	x <- x, y <- t.
	-- - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
	-- 		x <- t, y <- y.
	-- - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0, 
	-- 		x <- t, y <- x.
    if fx < ft then
        -- Case a
        y = t
        fy = ft
		dy = dt
	else
        -- Case c
        if dsign ~= 0 then
            y = x
            fy = fx
            dy = dx
		end
        -- Cases b and c
        x = t
        fx = ft
        dx = dt
    end

    -- Clip the new trial value in [tmin, tmax].
    if tmax < newt then newt = tmax end
    if newt < tmin then newt = tmin end

	-- Redefine the new trial value if it is close to the upper bound
	-- of the interval.
    if brackt ~= 0 and bound ~= 0 then
        mq = x + 0.66 * (y - x)
        if x < y then
            if mq < newt then newt = mq end
		else
            if newt < mq then newt = mq end
        end
    end

    -- Return the new trial value.
    t = newt
    return 0, x, fx, dx, y, fy, dy, t, brackt
end

local function line_search_morethuente(n, x, f, g, s, stp, xp, gp, wa, evaluate, param)
	local count = 0
	local uinfo = 0

    -- Check the input parameters for errors.
    if stp <= 0. then
		return "A logic error (negative line-search step) occurred.", f, stp
	end

    -- Compute the initial gradient in the search direction.
    local dginit = vecdot(g, s, n)

    -- Make sure that s points to a descent direction.
    if 0 < dginit then
        return lib.ERR_INCREASEGRADIENT, f, stp
	end

    -- Initialize local variables.
    local brackt = 0
    local stage1 = 1
    local finit = f
    local dgtest = param.ftol * dginit
    local width = param.max_step - param.min_step
    local prev_width = 2.0 * width

	-- The variables stx, fx, dgx contain the values of the step,
	-- function, and directional derivative at the best step.
	-- The variables sty, fy, dgy contain the value of the step,
	-- function, and derivative at the other endpoint of
	-- the interval of uncertainty.
	-- The variables stp, f, dg contain the values of the step,
	-- function, and derivative at the current step.
    local stx, sty = 0., 0.
    local fx, fy = finit, finit
    local dgx, dgy = dginit, dginit

	while true do
		-- Set the minimum and maximum steps to correspond to the
		-- present interval of uncertainty.
		local stmin, stmax
        if brackt ~= 0 then
            stmin = math.min(stx, sty)
            stmax = math.max(stx, sty)
		else
            stmin = stx
            stmax = stp + 4.0 * (stp - stx)
        end

        -- Clip the step in the range of [stpmin, stpmax].
        if stp < param.min_step then stp = param.min_step end
        if param.max_step < stp then stp = param.max_step end

		-- If an unusual termination is to occur then let
		-- stp be the lowest point obtained so far.
		if (brackt ~= 0 and ((stp <= stmin or stmax <= stp) or param.max_linesearch <= count + 1 or uinfo ~= 0)) or 
			(brackt ~= 0 and (stmax - stmin <= param.xtol * stmax)) then
            stp = stx
		end

		-- Compute the current value of x:
		--     x <- x + (*stp) * s.
        veccpy(x, xp, n)
        vecadd(x, s, stp, n)

        -- Evaluate the function and gradient values.
        f = evaluate(x, g, n, stp)
        local dg = vecdot(g, s, n)

		local ftest1 = finit + stp * dgtest
		count = count + 1

        -- Test for errors and convergence.
		if brackt ~= 0 and ((stp <= stmin or stmax <= stp) or uinfo ~= 0) then
			-- Rounding errors prevent further progress.
            return lib.ERR_ROUNDING_ERROR, f, stp
		end
        if stp == param.max_step and f <= ftest1 and dg <= dgtest then
            -- The step is the maximum value.
            return lib.ERR_MAXIMUMSTEP, f, stp
		end
        if stp == param.min_step and (ftest1 < f or dgtest <= dg) then
            -- The step is the minimum value.
            return lib.ERR_MINIMUMSTEP, f, stp
		end
        if brackt ~= 0 and (stmax - stmin) <= param.xtol * stmax then
            -- Relative width of the interval of uncertainty is at most xtol.
            return lib.ERR_WIDTHTOOSMALL, f, stp
		end
        if param.max_linesearch <= count then
            -- Maximum number of iteration.
            return lib.ERR_MAXIMUMLINESEARCH, f, stp
		end
        if f <= ftest1 and math.abs(dg) <= param.gtol * (-dginit) then
            -- The sufficient decrease condition and the directional derivative condition hold.
            return count, f, stp
		end

		-- In the first stage we seek a step for which the modified
		-- function has a nonpositive value and nonnegative derivative.
        if stage1 ~= 0 and f <= ftest1 and math.min(param.ftol, param.gtol) * dginit <= dg then
            stage1 = 0
		end

		-- A modified function is used to predict the step only if
		-- we have not obtained a step for which the modified
		-- function has a nonpositive function value and nonnegative
		-- derivative, and if a lower function value has been
		-- obtained but the decrease is not sufficient.
        if stage1 ~= 0 and ftest1 < f and f <= fx then
            -- Define the modified function and derivative values.
            local fm = f - stp * dgtest
            local fxm = fx - stx * dgtest
            local fym = fy - sty * dgtest
            local dgm = dg - dgtest
            local dgxm = dgx - dgtest
            local dgym = dgy - dgtest

			-- Call update_trial_interval() to update the interval of
			-- uncertainty and to compute the new step.
            uinfo, stx, fxm, dgxm, sty, fym, dgym, stp, brackt = update_trial_interval(
                stx, fxm, dgxm,
                sty, fym, dgym,
                stp, fm, dgm,
                stmin, stmax, brackt
            )

            -- Reset the function and gradient values for f.
            fx = fxm + stx * dgtest
            fy = fym + sty * dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
        else
			-- Call update_trial_interval() to update the interval of
			-- uncertainty and to compute the new step.
            uinfo, stx, fx, dgx, sty, fy, dgy, stp, brackt = update_trial_interval(
                stx, fx, dgx,
                sty, fy, dgy,
                stp, f, dg,
                stmin, stmax, brackt
			)
		end

		-- Force a sufficient decrease in the interval of uncertainty.
        if brackt ~= 0 then
            if 0.66 * prev_width <= math.abs(sty - stx) then
                stp = stx + 0.5 * (sty - stx)
			end
            prev_width = width
            width = math.abs(sty - stx)
        end
    end

	return lib.ERR_LOGICERROR, f, stp
end

---
-- Start a L-BFGS optimization.
--
--  @param  x           The array of variables. A client program can set
--                      default values for the optimization and receive the
--                      optimization result through this array.
--  @param  evaluate    The callback function to provide function and
--                      gradient evaluations given a current values of
--                      variables. A client program must implement a
--                      callback function compatible with \ref
--                      lbfgs_evaluate_t and pass the pointer to the
--                      callback function.
--  @param  progress    The callback function to receive the progress
--                      (the number of iterations, the current value of
--                      the objective function) of the minimization
--                      process. This argument can be set to \c nil or omit if
--                      a progress report is unnecessary.
--  @param  param       The table representing parameters for
--                      L-BFGS optimization. A client program can set this
--                      parameter to \c nil or omit to use the default parameters.
--  @return int         The status code. This function returns zero if the
--                      minimization process terminates without an error. A
--                      non-zero value indicates an error.
--  @return number      Same as x.
--  @param  number      The final value of the objective function for the variables.
local function lbfgs(x, evaluate, progress, param)
	if type(progress) == "table" then
		param = progress
		progress = nil
	else
		param = param or {}
	end
	local n = #x

	for k, v in pairs(defparam) do
		if param[k] == nil then param[k] = v end
	end
	local m = param.m

	local linesearch = line_search_morethuente

    -- Check the input parameters for errors.
	if n <= 0 then
		return lib.ERR_INVALID_N
	end
	if param.epsilon < 0. then
		return lib.ERR_INVALID_EPSILON
	end
	if param.past < 0 then
		return lib.ERR_INVALID_TESTPERIOD
	end
	if param.delta < 0. then
		return lib.ERR_INVALID_DELTA
	end
	if param.min_step < 0. then
		return lib.ERR_INVALID_MINSTEP
	end
	if param.max_step < param.min_step then
		return lib.ERR_INVALID_MAXSTEP
	end
	if param.ftol < 0. then
		return lib.ERR_INVALID_FTOL
	end
    if param.linesearch == lib.LINESEARCH_BACKTRACKING_WOLFE or
        param.linesearch == lib.LINESEARCH_BACKTRACKING_STRONG_WOLFE then
		if param.wolfe <= param.ftol or 1. <= param.wolfe then
			return lib.ERR_INVALID_WOLFE
		end
    end
	if param.gtol < 0. then
		return lib.ERR_INVALID_GTOL
	end
	if param.xtol < 0. then
		return lib.ERR_INVALID_XTOL
	end
	if param.max_linesearch <= 0 then
		return lib.ERR_INVALID_MAXLINESEARCH
	end
	if param.orthantwise_c < 0. then
		return lib.ERR_INVALID_ORTHANTWISE
	end
	if param.orthantwise_start < 1 or n < param.orthantwise_start then
		return lib.ERR_INVALID_ORTHANTWISE_START
	end
    if param.orthantwise_end < 1 then
        param.orthantwise_end = n;
	end
	if n < param.orthantwise_end then
		return lib.ERR_INVALID_ORTHANTWISE_END
	end
	if param.orthantwise_c ~= 0. then
		if param.linesearch == lib.LINESEARCH_BACKTRACKING then
			linesearch = line_search_backtracking_owlqn
		else
			-- Only the backtracking method is available.
			return lib.ERR_INVALID_LINESEARCH
		end
	else
		if param.linesearch == lib.LINESEARCH_MORETHUENTE then
			linesearch = line_search_morethuente
		elseif param.linesearch == lib.LINESEARCH_BACKTRACKING_ARMIJO or
			param.linesearch == lib.LINESEARCH_BACKTRACKING_WOLFE or
			param.linesearch == lib.LINESEARCH_BACKTRACKING_STRONG_WOLFE then
			linesearch = line_search_backtracking
		else
			return lib.ERR_INVALID_LINESEARCH
		end
	end

	local xp, g, gp, pg, d, w, lm = {}, {}, {}, {}, {}, {}, {}
	local pf = 0 < param.past and {} or nil
	for i = 1, m do
		lm[i] = {
			alpha = 0,
			ys = 0,
			s = {},
			y = {}
		}
	end

    -- Evaluate the function value and its gradient.
    local fx = evaluate(x, g, n, 0)
    if 0. ~= param.orthantwise_c then
        -- Compute the L1 norm of the variable and add it to the object value.
        local xnorm = owlqn_x1norm(x, param.orthantwise_start, param.orthantwise_end)
        fx = fx + xnorm * param.orthantwise_c
        owlqn_pseudo_gradient(
            pg, x, g, n,
            param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
            )
	end

	-- Store the initial value of the objective function.
	if pf then
		pf[0] = fx
	end

	-- Compute the direction;
	-- we assume the initial hessian matrix H_0 as the identity matrix.
    if param.orthantwise_c == 0. then
        vecncpy(d, g, n)
	else
        vecncpy(d, pg, n)
    end

	-- Make sure that the initial variables are not a minimizer.
	local xnorm = vec2norm(x, n)
	local gnorm = vec2norm(param.orthantwise_c == 0. and g or pg, n)
    if xnorm < 1.0 then xnorm = 1.0 end
	if gnorm / xnorm <= param.epsilon then
		return lib.ALREADY_MINIMIZED, x, fx
    end

    -- Compute the initial step:
    --     step = 1.0 / sqrt(vecdot(d, d, n))
    local step = vec2norminv(d, n)

    local k = 1
	local end_ = 0
	while true do
        -- Store the current position and gradient vectors.
        veccpy(xp, x, n)
        veccpy(gp, g, n)

		-- TODO
		local ls
        if param.orthantwise_c == 0. then
            ls, fx, step = linesearch(n, x, fx, g, d, step, xp, gp, w, evaluate, param)
		else
            ls, fx, step = linesearch(n, x, fx, g, d, step, xp, pg, w, evaluate, param)
            owlqn_pseudo_gradient(
                pg, x, g, n,
                param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
                )
		end
        if ls < 0 then
            -- Revert to the previous point.
            veccpy(x, xp, n)
            veccpy(g, gp, n)
			return ls, x, fx
		end

        -- Compute x and g norms.
		local xnorm = vec2norm(x, n);
		local gnorm = vec2norm(param.orthantwise_c == 0. and g or pg, n)

        -- Report the progress.
		if progress then
			local ret = progress(x, g, fx, xnorm, gnorm, step, n, k, ls)
			if ret and ret ~= 0 then 
				return ret, x, fx
			end
        end

		-- Convergence test.
		-- The criterion is given by the following formula:
		-- 	|g(x)| / \max(1, |x|) < \epsilon
        if xnorm < 1.0 then xnorm = 1.0 end
        if gnorm / xnorm <= param.epsilon then
            return lib.SUCCESS, x, fx
		end

		-- Test for stopping criterion.
		-- The criterion is given by the following formula:
		-- 	|(f(past_x) - f(x))| / f(x) < \delta
        if pf then
            -- We don't test the stopping criterion while k < past.
            if param.past <= k then
                -- Compute the relative improvement from the past.
                local rate = (pf[k % param.past] - fx) / fx

                -- The stopping criterion.
                if math.abs(rate) < param.delta then
					return lib.STOP, x, fx
				end
            end

            -- Store the current value of the objective function.
            pf[k % param.past] = fx
        end

        if param.max_iterations ~= 0 and param.max_iterations < k + 1 then
			-- Maximum number of iterations.
			return lib.ERR_MAXIMUMITERATION, x, fx
		end

		-- Update vectors s and y:
		--     s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
		--     y_{k+1} = g_{k+1} - g_{k}.
		local it = lm[end_ + 1]
        vecdiff(it.s, x, xp, n);
        vecdiff(it.y, g, gp, n);

		-- Compute scalars ys and yy:
		--     ys = y^t \cdot s = 1 / \rho.
		--     yy = y^t \cdot y.
		-- Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
        local ys = vecdot(it.y, it.s, n);
        local yy = vecdot(it.y, it.y, n);
        it.ys = ys;

		-- Recursive formula to compute dir = -(H \cdot g).
		--     This is described in page 779 of:
		--     Jorge Nocedal.
		--     Updating Quasi-Newton Matrices with Limited Storage.
		--     Mathematics of Computation, Vol. 35, No. 151,
		--     pp. 773--782, 1980.
        local bound = m <= k and m or k
        k = k + 1
        end_ = (end_ + 1) % m;

        -- Compute the steepest direction.
        if param.orthantwise_c == 0. then
            -- Compute the negative of gradients.
            vecncpy(d, g, n)
		else
            vecncpy(d, pg, n)
        end

		local j = end_
		for i = 1, bound do
            j = (j + m - 1) % m    -- if (--j == -1) j = m-1;
            local it = lm[j + 1]
            -- \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}.
            it.alpha = vecdot(it.s, d, n) / it.ys
            -- q_{i} = q_{i+1} - \alpha_{i} y_{i}.
            vecadd(d, it.y, -it.alpha, n)
        end

        vecscale(d, ys / yy, n)

        for i = 1, bound do
            local it = lm[j + 1]
            -- /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}.
            local beta = vecdot(it.y, d, n) / it.ys
            -- /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}.
            vecadd(d, it.s, it.alpha - beta, n)
            j = (j + 1) % m        -- if (++j == m) j = 0;
		end

		-- Constrain the search direction for orthant-wise updates.
        if param.orthantwise_c ~= 0. then
            for i = param.orthantwise_start, param.orthantwise_end do
                if d[i] * pg[i] >= 0 then
                    d[i] = 0
				end
            end
        end

		-- Now the search direction d is ready. We try step = 1 first.
        step = 1.0
    end
end

local strerrortb = {
	[lib.SUCCESS] = "Success: reached convergence (gtol).",
	[lib.STOP] = "Success: met stopping criteria (ftol).",
	[lib.ALREADY_MINIMIZED] = "The initial variables already minimize the objective function.",
	[lib.ERR_UNKNOWNERROR] = "Unknown error.",
	[lib.ERR_LOGICERROR] = "Logic error.",
	[lib.ERR_OUTOFMEMORY] = "Insufficient memory.",
	[lib.ERR_CANCELED] = "The minimization process has been canceled.",
	[lib.ERR_INVALID_N] = "Invalid number of variables specified.",
	[lib.ERR_INVALID_N_SSE] = "Invalid number of variables (for SSE) specified.",
	[lib.ERR_INVALID_X_SSE] = "The array x must be aligned to 16 (for SSE).",
	[lib.ERR_INVALID_EPSILON] = "Invalid parameter epsilon specified.",
	[lib.ERR_INVALID_TESTPERIOD] = "Invalid parameter past specified.",
	[lib.ERR_INVALID_DELTA] = "Invalid parameter delta specified.",
	[lib.ERR_INVALID_LINESEARCH] = "Invalid parameter linesearch specified.",
	[lib.ERR_INVALID_MINSTEP] = "Invalid parameter max_step specified.",
	[lib.ERR_INVALID_MAXSTEP] = "Invalid parameter max_step specified.",
	[lib.ERR_INVALID_FTOL] = "Invalid parameter ftol specified.",
	[lib.ERR_INVALID_WOLFE] = "Invalid parameter wolfe specified.",
	[lib.ERR_INVALID_GTOL] = "Invalid parameter gtol specified.",
	[lib.ERR_INVALID_XTOL] = "Invalid parameter xtol specified.",
	[lib.ERR_INVALID_MAXLINESEARCH] = "Invalid parameter max_linesearch specified.",
	[lib.ERR_INVALID_ORTHANTWISE] = "Invalid parameter orthantwise_c specified.",
	[lib.ERR_INVALID_ORTHANTWISE_START] = "Invalid parameter orthantwise_start specified.",
	[lib.ERR_INVALID_ORTHANTWISE_END] = "Invalid parameter orthantwise_end specified.",
	[lib.ERR_OUTOFINTERVAL] = "The line-search step went out of the interval of uncertainty.",
	[lib.ERR_INCORRECT_TMINMAX] = "A logic error occurred; alternatively, the interval of uncertainty" ..
		" became too small.",
	[lib.ERR_ROUNDING_ERROR] = "A rounding error occurred; alternatively, no line-search step" ..
		" satisfies the sufficient decrease and curvature conditions.",
	[lib.ERR_MINIMUMSTEP] = "The line-search step became smaller than min_step.",
	[lib.ERR_MAXIMUMSTEP] = "The line-search step became larger than max_step.",
	[lib.ERR_MAXIMUMLINESEARCH] = "The line-search routine reaches the maximum number of evaluations.",
	[lib.ERR_MAXIMUMITERATION] = "The algorithm routine reaches the maximum number of iterations.",
	[lib.ERR_WIDTHTOOSMALL] = "Relative width of the interval of uncertainty is at most" ..
		" xtol.",
	[lib.ERR_INCREASEGRADIENT] = "The current search direction increases the objective function value.",
}

local function strerror(err)
	return strerrortb[err] or "(unknown)"
end

lib.lbfgs = lbfgs
lib.strerror = strerror
return lib
