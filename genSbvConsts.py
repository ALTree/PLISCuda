import sys

if len(sys.argv) < 3:
    print("usage: python3 $genSbvConsts.py <sbc> filename.txt")
    sys.exit()

sbc = int(sys.argv[1])

with open(sys.argv[2], 'w') as f:
    for i in range(0, sbc):
        f.write(str(i) + ": 0\n")
    

