import sys

sbc = int(sys.argv[1])

with open(sys.argv[2], 'w') as f:
    for i in range(0, sbc):
        f.write(str(i) + ": 0\n")
    

