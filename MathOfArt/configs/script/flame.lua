_maxi = arg[1]
_width = arg[2]
_height = arg[3]
_scale = arg[4]

function hypot(x, y)
    return math.sqrt(x*x + y*y)
end

function gauss(mean, variance)
    return  math.sqrt(-2 * variance * math.log(math.random())) *
            math.cos(2 * math.pi * math.random()) + mean
end


function linear (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return x, y
end

function sinusoidal (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.sin(x), math.sin(y)
end

function spherical (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = 1.0 / (x * x + y * y)
    return r * x, r * y
end

function swirl (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = (x * x) + (y * y)
    return x * math.sin(r) - y * math.cos(r), x * math.cos(r) + y * math.sin(r)
end

function horseshoe (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = 1.0 / hypot(x, y)
    return r * (x - y) * (x + y), r * 2.0 * x * y
end

function polar (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.atan2(y, x) / pi, hypot(x, y) - 1.0
end

function handkerchief (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return r * math.sin(theta + r), r * math.cos(theta - r)
end

function heart (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return r * math.sin(theta * r), -r * math.cos(theta * r)
end

function disk (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y) * pi
    theta = math.atan2(y, x) / pi
    return theta * math.sin(r), theta * math.cos(r)
end

function spiral (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return (1.0 / r) * (math.cos(theta) + math.sin(r)), (1.0 / r) * (math.sin(theta) - math.cos(r))
end

function hyperbolic (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return math.sin(theta) / r, r * math.cos(theta)
end

function diamond (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return math.sin(theta) * math.cos(r), math.cos(theta) * math.sin(r)
end

function ex (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    i = math.sin(theta + r)
    i = i * i * i
    j = math.cos(theta - r)
    j = j * j * j
    return r * (i + j), r * (i - j)
end

function julia (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = math.sqrt(hypot(x, y))
    theta = math.atan2(y, x) * 0.5
    if math.random() > 0.5 then
        theta = theta + pi
    end
    return r * math.cos(theta), r * math.sin(theta)
end


function bent (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    if x >= 0.0 and y >= 0.0 then
        return x, y
    end
    if x < 0.0 and y >= 0.0 then
        return 2.0 *x, y
    end
    if x >= 0.0 and y < 0.0 then
        return x, y * 0.5
    end
    return 2.0 * x, y * 0.5
end

function waves (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return x + pa1 * math.sin(y / (pa2 * pa2)), y + pa3 * math.sin(x / (pa4 * pa4))
end

function fisheye (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = 2.0 / (1.0 + hypot(x, y))
    return r * y, r * x
end

function popcorn (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return x + c * math.sin(math.tan (3.0 * y)), y + f * math.sin(math.tan (3.0 * x))
end

function exponential (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.exp(x - 1.0) * math.cos(pi * y), math.exp(x - 1.0) * math.sin(pi * y)
end

function power (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    return math.pow(r, math.sin(theta)) * math.cos(theta), math.pow(r, math.sin(theta)) * math.sin(theta)
end

function cosine (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.cos(pi * x) * math.cosh (y), -math.sin(pi * x) * math.sinh (y)
end

function rings (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    p = pa2 * pa2
    prefix =((r + p) % (2.0 * p)) - p + (r * (1.0 - p))
    return prefix * math.cos(theta), prefix * math.sin(theta)
end

function fan (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    t = pi * c * c / 2 + sys.float_info.epsilon
    if (theta % (t * 2)) > t then
        return r * math.cos(theta - t), r * math.sin(theta - t)
    else
        return r * math.cos(theta + t), r * math.sin(theta + t)
    end
end

function eyefish (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = 2.0 / (1.0 + hypot(x, y))
    return r * x, r * y
end

function bubble (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = 4 + x * x + y * y
    return (4.0 * x) / r, (4.0 * y) / r
end

function cylinder (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.sin(x), y
end

function tangent (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.sin(x) / math.cos(y), math.tan(y)
end

function cross (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = math.sqrt (1.0 / ((x * x - y * y) * (x * x - y * y)))
    return x * r, y * r
end
function collatz (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return 0.25 * (1.0 + 4.0 * x - (1.0 + 2.0 * x) * math.cos(pi * x)), 0.25 * (1.0 + 4.0 * y - (1.0 + 2.0 * y) * math.cos(pi * y))
end
function mobius (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    t = (pa3 * x + pa4) * (pa3 * x + pa4) + pa3 * y * pa3 * y
    return ((pa1 * x + pa2) * (pa3 * x + pa4) + pa1 * pa3 * y * y) / t, (pa1 * y * (pa3 * x + pa4) - pa3 * y * (pa1 * x + pa2)) / t
end
function blob (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    theta = math.atan2(y, x)
    newx = r * (pa2 + 0.5 * (pa1 - pa2) * (math.sin(pa3 * theta) + 1)) * math.cos(theta)
    newy = r * (pa2 + 0.5 * (pa1 - pa2) * (math.sin(pa3 * theta) + 1)) * math.sin(theta)
    return nwx, newy
end
function noise (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    theta = math.random()
    r = math.random()
    return theta * x * math.cos(2 * pi * r), theta * y * math.sin(2 * pi * r)
end
function blur (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    theta = math.random()
    r = math.random()
    return theta * math.cos(2 * pi * r), theta * math.sin(2 * pi * r)
end

function square (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return math.random() - 0.5, math.random() - 0.5
end

function notBrokenWaves (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    return x + b * math.sin(y / math.pow(c, 2.0)), y + e * math.sin(x / math.pow(f, 2.0))
end

function juliaN (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4) --pa1 = power, pa2 = dist
    r = hypot(x, y)
    p3 = int(math.abs(pa1) * math.random())
    t = (math.atan2(x, y) + 2 * pi * p3)/pa1
    return math.pow(r, (pa2/pa1)) * math.cos(t), math.pow(r, (pa2/pa1)) * math.sin(t)
end

function juliaScope (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    r = hypot(x, y)
    p3 = math.floor(math.abs(pa1) * math.random())
    val = -1
    if math.random() > 0.5 then
        val = 1
    end
    t = ((val)*math.atan2(x, y) + 2 * math.pi * p3)/pa1
    return math.pow(r, (pa2/pa1)) * math.cos(t), math.pow(r, (pa2/pa1)) * math.sin(t)
end

function curl (x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4)
    t1 = 1 + (pa1*x) + (pa2 * (math.pow(x,2) - math.pow(y,2)))
    t2 = (pa1 * y) + (2 * pa2 * x * y)
    return (1/(math.pow(t1,2) + math.pow(t2, 2))) * (x*t1 + y*t2), (1/(math.pow(t1, 2) + math.pow(t2, 2))) * (y*t1 - x*t2)
end

function gaussian(x,y,a,b,c,d,e,f,pa1,pa2,pa3,pa4) --pa1 = mu, pa2 = sigma
    t = gauss(pa1, pa2)
    r = 2* pi * math.random()
    return t * math.cos(r), t * math.sin(r)
end
-----------------------------------------------------------------------------------


function rand(hi, lo)
    return lo + ((hi - lo) * math.random())
end

function sum(arr, sz)
    res = 0
    io.write(sz)
    for i = 0, sz do
        res = res + arr[i]
    end
    io.write(res)
    return res
end

function itterate(curX, curY, id)
    weight = weights[id]
    a = coefs[id][0]
    b = coefs[id][1]
    c = coefs[id][2]
    d = coefs[id][3]
    e = coefs[id][4]
    f = coefs[id][5]
    pa1 = params[id][0]
    pa2 = params[id][1]
    pa3 = params[id][2]
    pa4 = params[id][3]
    x = a * curX + b * curY + c
    y = d * curX + e * curY + f
    x, y = functions[id](x, y, a, b, c, d, e, f, pa1, pa2, pa3, pa4)
    curX = weight * x
    curY = weight * y
    return curX, curY
end

xres = 1920
yres = 1080
sup = 1
samples = 20000
seed = 23123
xmin = -1.777
xmax = 1.777
ymin = -1
ymax = 1
itters = 5000
functions = {}
functions[0]=ex
functions[1]=bubble
functions[2]=popcorn
functions[3]=waves
functions[4]=curl

probs = {}
probs[0] = 2
probs[1]=1
probs[2]=0.5
probs[3]=1.0
probs[4]=1.2

coefs = {}
for i = 0,5 do
    coefs[i] = {}
    for j = 0, 6 do
        coefs[i][j] = 0
    end
end

coefs[0][0] = 1
coefs[0][1] = 0
coefs[0][2] = 0.1
coefs[0][3] = 0
coefs[0][4] = 1
coefs[0][5] = 0.1

coefs[1][0] = 1
coefs[1][1] = 0.1
coefs[1][2] = 0
coefs[1][3] = 0.1
coefs[1][4] = 1
coefs[1][5] = 0

coefs[2][0] = 0
coefs[2][1] = 0
coefs[2][2] = 0
coefs[2][3] = 0
coefs[2][4] = 0
coefs[2][5] = 0

coefs[3][0] = 0
coefs[3][1] = 0
coefs[3][2] = 0
coefs[3][3] = 0
coefs[3][4] = 0
coefs[3][5] = 0

coefs[4][0] = -0.695157 
coefs[4][1] =  0.00633793
coefs[4][2] = 0.0475732
coefs[4][3] = -0.644204
coefs[4][4] = -0.270785
coefs[4][5] = 0.892879


params = {}
for i = 0,5 do
    params[i] = {}
end
params[0][0] = 0
params[0][1] = 0
params[0][2] = 0
params[0][3] = 0

params[1][0] = 2
params[1][1] = 1
params[1][2] = 0
params[1][3] = 0

params[2][0] = 0
params[2][1] = 0
params[2][2] = 0
params[2][3] = 0

params[3][0] = 0.269271
params[3][1] = -0.287245
params[3][2] = 1.40241
params[3][3] = -0.169667

params[4][0] = 0
params[4][1] = 0.4
params[4][2] = 0
params[4][3] = 0

weights = {}
weights[0]=1.0
weights[1]=1.0
weights[2]=1.0
weights[3]=0.8
weights[4]=1.2

gamma = 1

grid = {}
for i = 0, yres do
    grid[i] = {}
    for j = 0, xres do
        grid[i][j] = 0
    end
end

if seed < 0 then
    math.randomseed(os.time())
else
    math.randomseed(seed)
end

probsSum = 5.7--sum(probs, 3)
for i = 0, samples do
    curX = rand(xmin, xmax)
    curY = rand(ymin, ymax)
    
    
    for step = 0, itters do
        r = rand(0.000001, probsSum)
        val = 0.0
        id = 0
        while (r > val) do
            val = val + probs[id]
            id = id+1
        end
        id = id-1
        a, curX, curY = pcall(itterate, curX, curY, id)
        if not a then
            io.write("break")
            break
        end
        if step >= 20 then
            if xmin < curX and curX < xmax and ymin < curY and curY < ymax then
                x = math.floor( (curX - xmin)/(xmax - xmin) * xres * sup )
                y = math.floor( (curY - ymin)/(ymax - ymin) * yres * sup )
                if 0 < x and x < xres * sup and 0 < y and y < yres * sup then
                    grid[y][x] = grid[y][x] + 1
                end
            end
        end
    end
end

io.write("[Lua] out")
file = io.open("attractorGrid.txt", "w")
for i = 0, yres do
	for j = 0, xres do
		file:write(grid[i][j])
		file:write(" ")
	end
	file:write("\n")
end
file:close()
