



from numpy import float64


if __name__ == "__main__":

    # Read a file
    numbers=list()
    numbersS=list()
    GCOEFC1=[[0] * 161 for i in range(161)]
    GCOEFS1=[[0] * 161 for i in range(161)]
    lineNo=0
    with open("fileIO.txt", "r") as in_file:
        in_line = in_file.readline()
        numbers.append(in_line)
        in_line = in_file.readline()
        numbers.append(in_line)
        while True:
            in_line = in_file.readline()
            if not in_line:
                break
            in_line = in_line[:-1]
            if in_line[5]=="C":
                GCOEFC1[int(in_line[14:17])][int(in_line[17:22])]=(float(in_line[25:40])*10**float(in_line[41:]))
    # 				numbers.append([in_line[0:10],in_line[14:17],in_line[17:22],float(float(in_line[25:40])*10**float(in_line[41:]))])
            else:
                GCOEFS1[int(in_line[14:17])][int(in_line[17:22])]=(float(in_line[25:40])*10**float(in_line[41:]))
    # 				numbersS.append([in_line[0:10],in_line[14:17],in_line[17:22],float(float(in_line[25:40])*10**float(in_line[41:]))])

    print(GCOEFC1[152][92])

    fo = open("foo.txt", "w")
    for i in range(161):
        for j in range(161):

            t=str(i)
            k=str(j)
            ok = (repr(GCOEFC1[(i)][(j)]))

            if (ok!="0"):
                fo.write("GCOEFC1"+" "+t+" "+k+" "+ok+" "+"\n")
    fo.close()



