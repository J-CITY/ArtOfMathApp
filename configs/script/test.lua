-- io.write("[Lua] These args were passed into the script from Cn")
-- for i=1,#arg do
--     print("      ", i, arg[i])
-- end 

zr = arg[1]
zi = arg[2]
cr = arg[3]
ci = arg[4]
maxi = arg[5]

-----------------------------------
local cmath =  {}
--module("cmath", package.seeall)

-- Logarithm

function log (x)
	local xx = tocomplex(x)
	if ((xx.re <= 0) and (xx.im == 0)) then
		if (xx.re == 0) then
			return new(-math.huge, 0)
		else
			return nil
		end
	else
		return new(math.log(xx:abs()), xx:arg())
	end
end

function log10 (x)
	local xx = tocomplex(x)
	return (log(xx) / log(10))
end

-- Exponent

function exp (x)
	local xx = tocomplex(x)
	if (xx.im == 0) then
		return new(math.exp(xx.re), 0)
	else
		return new(
			math.exp(xx.re) * math.cos(xx.im),
			math.exp(xx.re) * math.sin(xx.im)
		)
	end
end

function pow (x, y)
	local base = tocomplex(x)
	local exponent = tocomplex(y)
	if (exponent.im ~= 0) then
		return exp(exponent * log(base))
	elseif ((base.re < 0) and (base.im == 0) and
		(exponent.re * 2 == math.floor(exponent.re * 2))) then
		if (math.floor(exponent.re * 2) % 4 == 0) then
			return new((-base.re) ^ exponent.re, 0)
		elseif (math.floor(exponent.re * 2) % 4 == 1) then
			return new(0, (-base.re) ^ exponent.re)
		elseif (math.floor(exponent.re * 2) % 4 == 2) then
			return new(-((-base.re) ^ exponent.re), 0)
		elseif (math.floor(exponent.re * 2) % 4 == 3) then
			return new(0, -((-base.re) ^ exponent.re))
		end
	elseif ((base.re == 0) and (base.im ~= 0) and
		(exponent.re == math.floor(exponent.re))) then
		if (math.floor(exponent.re) % 4 == 0) then
			return new(base.im ^ exponent.re, 0)
		elseif (math.floor(exponent.re) % 4 == 1) then
			return new(0, base.im ^ exponent.re)
		elseif (math.floor(exponent.re) % 4 == 2) then
			return new(-(base.im ^ exponent.re), 0)
		elseif (math.floor(exponent.re) % 4 == 3) then
			return new(0, -(base.im ^ exponent.re))
		end
	elseif (base:arg() ~= 0) then
		return polar(
			base:abs() ^ exponent.re,
			base:arg() * exponent.re
		)
	else
		return new(base.re ^ exponent.re, 0)
	end
end

function sqrt (x)
	local xx = tocomplex(x)
	return pow(xx, 0.5)
end

-- Trigonometric

function sin (x)
	local xx = tocomplex(x)
	return new(
		math.sin(xx.re) * math.cosh(xx.im),
		math.cos(xx.re) * math.sinh(xx.im)
	)
end

function cos (x)
	local xx = tocomplex(x)
	return new(
		math.cos(xx.re) * math.cosh(xx.im),
		-(math.sin(xx.re)) * math.sinh(xx.im)
	)
end

function tan (x)
	local xx = tocomplex(x)
	return (sin(xx) / cos(xx))
end

function asin (x)
	local xx = tocomplex(x)
	return (-(i) * log(i * xx + sqrt(1 - xx ^ 2)))
end

function acos (x)
	local xx = tocomplex(x)
	return (math.pi / 2 - asin(xx))
end

function atan (x)
	local xx = tocomplex(x)
	return (0.5 * i * log((1 - i * xx) / (1 + i * xx)))
end

-- Hyperbolic

function sinh (x)
	local xx = tocomplex(x)
	return ((exp(x) - exp(-x)) / 2)
end

function cosh (x)
	local xx = tocomplex(x)
	return ((exp(x) + exp(-x)) / 2)
end

function tanh (x)
	local xx = tocomplex(x)
	return (sinh(xx) / cosh(xx))
end

function asinh (x)
	local xx = tocomplex(x)
	return log(xx + sqrt(xx ^ 2 + 1))
end

function acosh (x)
	local xx = tocomplex(x)
	return log(xx + sqrt(xx + 1) * sqrt(xx - 1))
end

function atanh (x)
	local xx = tocomplex(x)
	return (0.5 * log((1 + xx) / (1 - xx)))
end

-- Complex number specific functions ()

function conj (x)
	local xx = tocomplex(x)
	return (xx:conj())
end
function abs (x)
	local xx = tocomplex(x)
	return (xx:abs())
end
function arg (x)
	local xx = tocomplex(x)
	return (xx:arg())
end



-- Module "complex" ... prototype for the complex number objects
local complex =  {}
--module("complex", package.seeall)

-- Assumes number as a complex number
function tocomplex (x)
	if (type(x) == "number") then return new(x, 0) else return x end
end

-- x:conj() ... Complex conjugate
function conj (x)
	return new(x.re, -(x.im))
end

-- x:abs() ... Norm
function abs (x)
	return math.sqrt(x.re ^ 2 + x.im ^ 2)
end

-- x:arg() ... Argument
function arg (x)
	return math.atan2(x.im, x.re)
end

-- Something like operator overloading
complex_meta = {
	-- Addition .. (a+bi)+(c+di) == (a+c)+(b+d)i
	__add = function (x, y)
		local xx = tocomplex(x); local yy = tocomplex(y)
		return new(xx.re + yy.re, xx.im + yy.im)
	end,

	-- Subtraction .. (a+bi)-(c+di) == (a+c)-(b+d)i
	__sub = function (x, y)
		local xx = tocomplex(x); local yy = tocomplex(y)
		return new(xx.re - yy.re, xx.im - yy.im)
	end,

	-- Unary minus .. -(a+bi) == -a-bi
	__unm = function (x)
		local xx = tocomplex(x)
		return new(-xx.re, -xx.im)
	end,

	-- Multiplication .. (a+bi)*(c+di) == (ac-bd)+(ad+bc)i
	__mul = function (x, y)
		local xx = tocomplex(x); local yy = tocomplex(y)
		return new(
			xx.re * yy.re - xx.im * yy.im,
			xx.re * yy.im + xx.im * yy.re
		)
	end,

	-- Division .. (a+bi)/(c+di) == ((ac+bd)+(bc-ad)i)/(c^2+d^2)
	__div = function (x, y)
		local xx = tocomplex(x); local yy = tocomplex(y)
		return new(
			(xx.re * yy.re + xx.im * yy.im) / (yy.re * yy.re + yy.im * yy.im),
			(xx.im * yy.re - xx.re * yy.im) / (yy.re * yy.re + yy.im * yy.im)
		)
	end,

	-- yth power of x
	__pow = function (x, y)
		local xx = tocomplex(x); local yy = tocomplex(y)
		return cmath.pow(xx, yy)
	end,

	-- mod and concat are unsupported
	__mod = function (x, y) error ("unsupported operator") end,
	__concat = function (x, y) error ("unsupported operator") end,

	-- Equality
	__eq = function (x, y)
		if ((x.re == y.re) and (x.im == y.im)) then
			return true
		else
			return false
		end
	end,

	-- Sorry, complex numbers are not comparable...
	__lt = function (x, y) error ("complex numbers are not comparable") end,
	__le = function (x, y) error ("complex numbers are not comparable") end,

	-- tostring()
	__tostring = function (x)
		if (x.re == 0) then
			if (x.im == 0) then
				return "0"
			else
				return "" .. x.im .. "i"
			end
		else
			if (x.im > 0) then
				return "" .. x.re .. "+" .. x.im .. "i"
			elseif (x.im < 0) then
				return "" .. x.re .. x.im .. "i"
			else
				return "" .. x.re
			end
		end
	end,

	-- tonumber() ... works only if it is actually a real number
	__tonumber = function(x)
		if (x.im == 0) then
			return x.re
		else
			return nil
		end
		
	end,
}

-- Constructor
function new (r, i)
	local cn = {re = r, im = i}
	setmetatable(cn, complex_meta)
	cn.conj = conj
	cn.abs = abs
	cn.arg = arg
	return cn
end

-- Another constructor ... by polar form
function polar (r, theta)
	return new(r * math.cos(theta), r * math.sin(theta))
end
----------------------------------
-----------------------------------
z = new(zr, zi)
c = new(cr, ci)

i = 0
while ((z.re*z.re + z.im*z.im < 4) and i < maxi) do
    -- WRITE CODE HERE
    -- z = ch(z)
    -- z = ctg(z)
    -- z = sin(z) * cos(z)
    -- z = sh(z)
    -- z = tg(z)
    -- z = z * cos(z)
    -- z = z * sin(z)
    z = z*z + c
    c = c/2 + z
    -- END YOUR CODE
    i = i + 1
end

[[ Mandelbrot (for cloud)
    z = <...>
    c = <...>
]]

[[ Julia 
i = 0
r = 0.7 -- same as in config.toml
while (i < maxi) do
    if (r > 0 and (math.sqrt(z.re*z.re + z.im*z.im) < r)) then break end
    z = <...>
    c = <...>
    i = i + 1
end
]]

[[ Julia (for cloud)
    z = <...>
    c = <...>
]]

newz = z
newc = c

return newz.re, newz.im, newc.re, newc.im, i