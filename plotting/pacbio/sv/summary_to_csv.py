import sys
with open(sys.argv[1], 'r') as f:
    hdr = []
    value = []
    for line in f:
        line = line.split(':')
        if len(line) > 1:
            l0 = line[0]
            l1 = line[1].rstrip()
            hdr.append(l0[l0.find('"')+1 : -1])
            value.append(l1[: l1.find(',')])
            # print(l0[l0.find('"')+1 : -1], l1[: l1.find(',')], sep=',')
    print(','.join(hdr))
    print(','.join(value))
