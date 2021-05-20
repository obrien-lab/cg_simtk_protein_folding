def ReadLastColumn(filename):
    QList = []
    results = []
    n = 0
    f = open(filename, "r")
    for line in f:
        if n == 0:
            n = n + 1
            continue
        sline = line.strip().split()
        QList.append(float(sline[-1]))
    f.close()
    results.append(QList)
    return results