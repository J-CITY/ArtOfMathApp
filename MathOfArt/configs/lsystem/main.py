f = open('all.toml', 'r')
fout = None

for line in f:
    if line.isspace():
        continue
    if line[0] == '#':
        if fout is not None:
            fout.close()
        name = line[1:].strip()
        name = name.replace(' ', '_')
        fout = open(name+'.toml', 'w')
        fout.write(line)
        fout.write("iterations = 2\n")
        fout.write("distance = 10\n")
        fout.write("line_width = 3\n")
        continue
    fout.write(line)



