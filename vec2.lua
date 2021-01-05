local lib = {}

function lib.add(a, b)
	return lib.new(a[1] + b[1], a[2] + b[2])
end

function lib.sub(a, b)
	return lib.new(a[1] - b[1], a[2] - b[2])
end

function lib.mul(a, b)
	return lib.new(a[1] * b, a[2] * b)
end

function lib.div(a, b)
	return lib.new(a[1] / b, a[2] / b)
end

function lib.dot(a, b)
	return a[1] * b[1] + a[2] * b[2]
end

function lib.len2(v)
	return lib.dot(v, v)
end

function lib.len(v)
	return math.sqrt(lib.len2(v))
end

function lib.angle(v)
	return math.atan2(v[2], v[1])
end

function lib.rot(a, b)
	local sinb, cosb = math.sin(b), math.cos(b)
	return lib.new(a[1] * cosb - a[2] * sinb, a[1] * sinb + a[2] * cosb)
end

function lib.eq(a, b)
	return a[1] == b[1] and a[2] == b[2]
end

function lib.tostring(v)
	return string.format("vec2 {%f, %f}", v[1], v[2])
end

function lib.toVec3(v)
	return { v[1], 0.0, v[2] }
end

local mt = {
	__index = lib,
	__add = lib.add,
	__sub = lib.sub,
	__mul = lib.mul,
	__div = lib.div,
	__eq = lib.eq,
	__tostring = lib.tostring,
}

function lib.new(x, y)
	local v = type(x) == "table" and { unpack(x) } or { x, y }
	setmetatable(v, mt)
	return v
end

function lib.fromVec3(t)
	return lib.new(t[1], t[3])
end

return lib
