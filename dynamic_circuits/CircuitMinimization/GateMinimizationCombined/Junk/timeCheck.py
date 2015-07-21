import time
def OR1():
    a = "00001111"
    b = "00110011"
    answer = ""
    for i in range(8):
        if (a[i]=="1" or b[i]=="1"):
            answer += "1"
        elif (a[i]=="0"  and b[i]=="0"):
            answer += "0"
    return answer
def OR2():
    a = "00001111"
    b = "00110011"
    answer = ""
    for i in range(8):
        answer += str(int(a[i]) or int(b[i]))
    return answer
def OR3():
    a = "00001111"
    b = "00110011"
    answer = ""
    for i in range(8):
        if (a[i]=="1" or b[i]=="1"):
            answer += "1"
        else:
            answer += "0"
    return answer
def OR4():
    a = "00001111"
    b = "00110011"
    answer = ""
    for i in range(8):
        if (a[i]=="0"  and b[i]=="0"):
            answer += "0"
        else:
            answer += "1"
    return answer
def compare(num):
    startTime = time.time()
    for i in range(num):
        x = OR1()
    endTime = time.time()
    print "Version 1 took",str(endTime-startTime),"seconds."
    startTime = time.time()
    for i in range(num):
        x = OR2()
    endTime = time.time()
    print "Version 2 took",str(endTime-startTime),"seconds."
    startTime = time.time()
    for i in range(num):
        x = OR3()
    endTime = time.time()
    print "Version 3 took",str(endTime-startTime),"seconds."
    startTime = time.time()
    for i in range(num):
        x = OR4()
    endTime = time.time()
    print "Version 4 took",str(endTime-startTime),"seconds."
