import random
import gmpy2
import hashlib
FK = [0xa3b1bac6, 0x56aa3350, 0x677d9197, 0xb27022dc]
CK = [0x00070e15, 0x1c232a31, 0x383f464d, 0x545b6269, 0x70777e85, 0x8c939aa1, 0xa8afb6bd, 0xc4cbd2d9,
      0xe0e7eef5, 0xfc030a11, 0x181f262d, 0x343b4249, 0x50575e65, 0x6c737a81, 0x888f969d, 0xa4abb2b9,
      0xc0c7ced5, 0xdce3eaf1, 0xf8ff060d, 0x141b2229, 0x30373e45, 0x4c535a61, 0x686f767d, 0x848b9299,
      0xa0a7aeb5, 0xbcc3cad1, 0xd8dfe6ed, 0xf4fb0209, 0x10171e25, 0x2c333a41, 0x484f565d, 0x646b7279]
Sbox = [[0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05],
        [0x2b, 0x67, 0x9a, 0x76, 0x2a, 0xbe, 0x04, 0xc3, 0xaa, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99],
        [0x9c, 0x42, 0x50, 0xf4, 0x91, 0xef, 0x98, 0x7a, 0x33, 0x54, 0x0b, 0x43, 0xed, 0xcf, 0xac, 0x62],
        [0xe4, 0xb3, 0x1c, 0xa9, 0xc9, 0x08, 0xe8, 0x95, 0x80, 0xdf, 0x94, 0xfa, 0x75, 0x8f, 0x3f, 0xa6],
        [0x47, 0x07, 0xa7, 0xfc, 0xf3, 0x73, 0x17, 0xba, 0x83, 0x59, 0x3c, 0x19, 0xe6, 0x85, 0x4f, 0xa8],
        [0x68, 0x6b, 0x81, 0xb2, 0x71, 0x64, 0xda, 0x8b, 0xf8, 0xeb, 0x0f, 0x4b, 0x70, 0x56, 0x9d, 0x35],
        [0x1e, 0x24, 0x0e, 0x5e, 0x63, 0x58, 0xd1, 0xa2, 0x25, 0x22, 0x7c, 0x3b, 0x01, 0x21, 0x78, 0x87],
        [0xd4, 0x00, 0x46, 0x57, 0x9f, 0xd3, 0x27, 0x52, 0x4c, 0x36, 0x02, 0xe7, 0xa0, 0xc4, 0xc8, 0x9e],
        [0xea, 0xbf, 0x8a, 0xd2, 0x40, 0xc7, 0x38, 0xb5, 0xa3, 0xf7, 0xf2, 0xce, 0xf9, 0x61, 0x15, 0xa1],
        [0xe0, 0xae, 0x5d, 0xa4, 0x9b, 0x34, 0x1a, 0x55, 0xad, 0x93, 0x32, 0x30, 0xf5, 0x8c, 0xb1, 0xe3],
        [0x1d, 0xf6, 0xe2, 0x2e, 0x82, 0x66, 0xca, 0x60, 0xc0, 0x29, 0x23, 0xab, 0x0d, 0x53, 0x4e, 0x6f],
        [0xd5, 0xdb, 0x37, 0x45, 0xde, 0xfd, 0x8e, 0x2f, 0x03, 0xff, 0x6a, 0x72, 0x6d, 0x6c, 0x5b, 0x51],
        [0x8d, 0x1b, 0xaf, 0x92, 0xbb, 0xdd, 0xbc, 0x7f, 0x11, 0xd9, 0x5c, 0x41, 0x1f, 0x10, 0x5a, 0xd8],
        [0x0a, 0xc1, 0x31, 0x88, 0xa5, 0xcd, 0x7b, 0xbd, 0x2d, 0x74, 0xd0, 0x12, 0xb8, 0xe5, 0xb4, 0xb0],
        [0x89, 0x69, 0x97, 0x4a, 0x0c, 0x96, 0x77, 0x7e, 0x65, 0xb9, 0xf1, 0x09, 0xc5, 0x6e, 0xc6, 0x84],
        [0x18, 0xf0, 0x7d, 0xec, 0x3a, 0xdc, 0x4d, 0x20, 0x79, 0xee, 0x5f, 0x3e, 0xd7, 0xcb, 0x39, 0x48]]
def egcd(a, b):
    #a*xi + b*yi = ri
    if b == 0:
        return (1, 0, a)
    #a*x1 + b*y1 = a
    x1 = 1
    y1 = 0
    #a*x2 + b*y2 = b
    x2 = 0
    y2 = 1
    while b != 0:
        q =int(a // b)
        #ri = r(i-2) % r(i-1)
        r = a % b
        a = b
        b = r
        #xi = x(i-2) - q*x(i-1)
        x = x1 - q*x2
        x1 = x2
        x2 = x
        #yi = y(i-2) - q*y(i-1)
        y = y1 - q*y2
        y1 = y2
        y2 = y
    if a < 0:
        a = -1 * a
        x1 = -1 * x1
        y1 = -1 * y1
    while x1 < 0:
        if(x2>0):
            x1=x2 +x1
            y1=y2 +y1
        else:
            x1 = x1 - x2
            y1 = y1 -y2
    while (x1-abs(x2)) > 0:
        if(x2<0):
            x1=x2 +x1
            y1=y2 +y1
        else:
            x1 = x1 - x2
            y1 = y1 -y2
    return(x1, y1, a)
def pow(b, e, m):
    result = 1
    while e != 0:
        if (e&1) == 1:
            # ei = 1, then mul
            result = (result * b) % m
        e >>= 1
        # b, b^2, b^4, b^8, ... , b^(2^n)
        b = (b*b) % m
    return result
def long2str(a):
    return int.to_bytes(a,len(bin(a)[2:])//8+1,byteorder='big')
def str2long(a):
    return int.from_bytes(a, byteorder='big')
def inverse(a,p):
    return egcd(a,p)[0]
def cht(a,m):
    M=1
    for k in m:
        M=M*k
    x=0
    for i in range(len(a)):
        t=egcd((M//m[i]),m[i])[0]
        x=(x+a[i]*t*(M//m[i]))%M
    return x
def check(num):# miller_rabin
    y = num - 1
    r = 0
    if num < 2: return False
    if num < 2: return False
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
                    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443,
                    449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
                    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
                    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
                    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983,
                    991, 997]
    if num in small_primes: return True  # 如果是小素数,那么直接返回true
    # 如果大数是这些小素数的倍数,那么就是合数,返回false
    for prime in small_primes:
        if num % prime == 0: return False

    while y % 2 == 0:
        y = y // 2
        r += 1

    for _ in range(10):
        a = random.randint(2, num - 1)
        while egcd(a, num)[2] != 1:
            a = random.randint(2, num - 1)
        x = pow(a, y, num)
        if x == 1 or x == num - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, num)
            if x == 1:
                return False
            if x == num - 1:
                break
        if x != num - 1:
            return False
    return True
def getPrime(n):
    i = random.randint(2 ** (n-1), 2 ** n)
    while (check(i) == 0):
        i = random.randint(2 ** (n-1), 2 ** n)
    return i
def left_move(input,n):
    list='{:032b}'.format(input)
    return int(list[n:len(list)]+list[0:n],2)
def SBox(input):
    small=[0 for i in range(0,4)]
    out=0
    for i in range(0,4):
        small[3-i]=(input >>(i*8))&0b11111111
    for i in range(4):
        l = small[i] & 0b1111
        h = (small[i] >> 4) & 0b1111
        out <<= 8
        out += Sbox[h][l]
    return out
def key_create(k):
    MK = [0 for i in range(0, 4)]
    K = [0 for i in range(0, 36)]
    for i in range(0, 4):
        MK[3-i]=(k>>(i*32))&0b11111111111111111111111111111111
    for i in range(0, 4):
        K[i] = MK[i] ^ FK[i]
    for i in range(0, 32):
        t = K[i + 1] ^ K[i + 2] ^ K[i + 3]
        t = t ^ CK[i]
        t = SBox(t)
        state = t ^ left_move(t, 13) ^ left_move(t, 23)
        K[i + 4] = state ^ K[i]
    return K
def sm4_encrypt(s,K):
    X = [0 for i in range(0, 36)]
    for i in range(0, 4):
        X[3 - i] = (s >> (i * 32)) & 0b11111111111111111111111111111111
    for i in range(0, 32):
        t = X[i + 1] ^ X[i + 2] ^ X[i + 3]
        t = t ^ K[i + 4]
        t = SBox(t)
        state = t ^ left_move(t, 2) ^ left_move(t, 10) ^ left_move(t, 18) ^ left_move(t, 24)
        X[i + 4] = state ^ X[i]
        # print(hex(X[i+4]))
    out=0
    for i in range(0,4):
        out=(out<<32)+X[35-i]
    return out
def sm4_decrypt(s,K):
    X = [0 for i in range(0, 36)]
    for i in range(0, 4):
        X[3 - i] = (s >> (i * 32)) & 0b11111111111111111111111111111111
    for i in range(0, 32):
        t = X[i + 1] ^ X[i + 2] ^ X[i + 3]
        t = t ^ K[31-i + 4]
        t = SBox(t)
        state = t ^ left_move(t, 2) ^ left_move(t, 10) ^ left_move(t, 18) ^ left_move(t, 24)
        X[i + 4] = state ^ X[i]
        # print(hex(X[i+4]))
    out=0
    for i in range(0,4):
        out=(out<<32)+X[35-i]
    return out
def SM4(s,op,K):
    if op==1:
        return sm4_encrypt(s,K)
    else:
        return sm4_decrypt(s,K)
def outprint(t):
    tmp=[0 for i in range(16)]
    out=""
    for i in range(16):
        tmp[15-i]=(t>>8*i)&0xff
    for i in tmp:
        out=out+'{:02x}'.format(i)
    return out
def CBC(k,M,iv,op):
    out=""
    P = []
    K = key_create(k)
    data = []
    for i in range(len(M)):
        # data[i] = int(str(tmp[i]), 16)
        data.append(M[i])
    sum = 0
    for i in range(len(data)):
        if i % 16 == 0 and i != 0:
            P.append(sum)
            sum = 0
        sum = (sum << 8) + data[i]
    if i % 16 == 15 and op == 1:
        P.append(sum)
    if i % 16 == 15 and op == 0:
        P.append(sum)
    if op == 1:
        for k in range(len(data) % 16, 16):
            sum = (sum << 8) + 16 - len(data) % 16
        P.append(sum)
    c = iv
    if op == 1:
        for i in range(len(P)):
            c = SM4(c ^ P[i], op, K)
            out=out+outprint(c)
    else:
        for i in range(len(P) - 1):
            p = c ^ SM4(P[i], op, K)
            c = P[i]
            out=out+outprint(p)
        ans = c ^ SM4(P[i + 1], op, K)
        temp = ans % 16
        if (temp != 0):
            ans = ans >> (8 * temp)
            tmp = [0 for i in range(16 - temp)]
            for i in range(16 - temp):
                tmp[15 - temp - i] = (ans >> 8 * i) & 0xff
            for i in tmp:
                out=out+ '{:02x}'.format(i)
            out=out+''
    return long2str(int(out,16))


def CTR(k,M,iv,op):
    t2 = [0 for i in range(0, 16)]
    t = []
    out = ""
    P = []
    K = key_create(k)
    data = []
    i = 0
    for i in range(len(M)):
        data.append(M[i])
    for i in range(iv, iv + len(data) // 16 + 1):
        q = SM4(i, 1, K)
        for j in range(0, 16):
            t2[j] = q & 0xff
            q >>= 2 * 4
        for j in range(0, 16):
            t.append(t2[15 - j])
    for i in range(len(data)):
        if i % 16 == 0:
            out=out+""
        out=out+ '{:02x}'.format(data[i] ^ t[i])
    return long2str(int(out,16))
def bytes_read(input_file):
    with open(input_file, 'rb') as f:
        f.seek(0, 2)  # 读取文件的总长度，seek(0,2)移到文件末尾，tell()指出当前位置，并且用seek(0)重新回到起点
        size = f.tell()
        f.seek(0)
        bytes_list = [0] * size  # 创建一个长度为size的列表，存放读入的字节

        i = 0
        while i < size:
            bytes_list[i] = f.read(1)  # 每次读取一个符号
            i += 1
    return bytes_list

def CBC_file_encrypt(inputfile_path,outputfile_path,k,iv):
    f=open(outputfile_path,'wb')
    bytes_list = bytes_read(inputfile_path)
    data=b""
    for q in bytes_list:
        data=data+q
    out=CTR(k,data,iv,1)
    f.write(out)
    f.close
def CBC_file_decrypt(inputfile_path,outputfile_path,k,iv):
    f=open(outputfile_path,'wb')
    bytes_list = bytes_read(inputfile_path)
    data = b""
    for q in bytes_list:
        data = data + q
    out=CTR(k,data,iv,0)
    f.write(out)
    f.close
def CTR_file_encrypt(inputfile_path,outputfile_path,k,iv):
    f=open(outputfile_path,'wb')
    bytes_list = bytes_read(inputfile_path)
    data=b""
    for q in bytes_list:
        data=data+q
    # print(data)
    out=CTR(k,data,iv,1)
    f.write(out)
    f.close
def CTR_file_decrypt(inputfile_path,outputfile_path,k,iv):
    f=open(outputfile_path,'wb')
    bytes_list = bytes_read(inputfile_path)
    data = b""
    for q in bytes_list:
        data = data + q
    out=CTR(k,data,iv,0)
    f.write(out)
    f.close
def RSA_encrypt(m,e,N):
    c=pow(m,e,N)
    return c
def RSA_decrypt(c,e,p,q):
    (x,y,z)=egcd((p-1)*(q-1),e)
  #  print(x,y,z)
    d=int(y)
    if(d<0):
        d=d+(p-1)*(q-1)
    # m1=((c%p)**(d%(p-1)))%p
    m1=pow(c%p,d%(p-1),p)
    #m2 = ((c % q) ** (d % (q - 1))) % q
    m2 = pow(c % q, d % (q - 1), q)
    A=egcd(q,p)[0]
    B=egcd(p,q)[0]
    m=(m1*A*q+m2*B*p)%(p*q)
    return m
def pad(ming):
    binText = str_2_bin(ming)
    # print(len(binText)-2**64)
    initLen = binLen = len(binText)
    binText += '1'
    binLen += 1
    while binLen % 512 != 448:
        binText += '0'
        binLen += 1
    temp = bin(initLen)[2:]
    lentemp = len(temp)
    while lentemp != 64:
        temp = '0' + temp
        lentemp += 1
    binText = binText + temp
    # print(binText)
    return binText


def padhex(ming):
    binText = hex_2_bin(ming)
    initLen = binLen = len(binText)
    binText += '1'
    binLen += 1
    while binLen % 512 != 448:
        binText += '0'
        binLen += 1
    temp = bin(initLen)[2:]
    lentemp = len(temp)
    while lentemp != 64:
        temp = '0' + temp
        lentemp += 1
    binText = binText + temp
    # print(binText)
    return binText


def hex_2_bin(strr):
    result = ''
    for i in strr:
        x = bin(int('0x'+i, 16))[2:]
        while len(x) < 4:
            x = '0' + x
        result += x
    return result


def str_2_bin(strr):
    result = ''
    for i in strr:
        x = bytes(i, encoding="utf-8")
        for j in x:
            y = bin(j)[2:]
            while len(y) % 8 != 0:
                y = '0' + y
            result += y
    return result


def Hash1(ming):
    h0 = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0]
    EM = padhex(ming)
    # print(len(EM))
    fenzu = len(EM) // 512
    M = []
    h = h0
    for i in range(fenzu):
        M.append(EM[512 * i:512 * (i + 1)])
    for j in range(len(M)):
        W = extend(M[j])
        h = Round(W, h)
    result = 0
    for k in range(len(h)):
        result = result << 32 | h[k]
    # print(hex(result))
    return result


def Hash(ming):
    h0 = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0]
    EM = pad(ming)
    # print(len(EM))
    fenzu = len(EM) // 512
    M = []
    h = h0
    for i in range(fenzu):
        M.append(EM[512 * i:512 * (i + 1)])
    for j in range(len(M)):
        W = extend(M[j])
        h = Round(W, h)
    result = 0
    for k in range(len(h)):
        result = result << 32 | h[k]
    return result


def extend(M):
    W = []
    for i in range(16):
        W.append(int(M[32 * i:32 * (i + 1)], 2))
    '''for i in range(16):
        print(bin(W[i]))'''
    for j in range(16, 80):
        temp = W[j - 3] ^ W[j - 8] ^ W[j - 14] ^ W[j - 16]
        temp = ((temp << 1) & 0b11111111111111111111111111111110) | temp >> 31
        # print(bin(temp))
        W.append(temp)
    return W


def Round(W, h0):
    a = h0[0]
    b = h0[1]
    c = h0[2]
    d = h0[3]
    e = h0[4]
    h = [a, b, c, d, e]
    for i in range(80):
        temp1 = ((a << 5) & 0b11111111111111111111111111100000) | a >> 27
        temp2 = f(b, c, d, i)
        temp3 = e
        temp4 = W[i]
        temp5 = K(i)
        temp = (temp1 + temp2 + temp3 + temp4 + temp5) % 2 ** 32
        e = d
        d = c
        c = ((b << 30) & 0b11000000000000000000000000000000) | b >> 2
        b = a
        a = temp
    '''print(bin(a))
    print(bin(b))
    print(bin(c))
    print(bin(d))
    print(bin(e))'''
    h[0] = (h[0] + a) % 2 ** 32
    h[1] = (h[1] + b) % 2 ** 32
    h[2] = (h[2] + c) % 2 ** 32
    h[3] = (h[3] + d) % 2 ** 32
    h[4] = (h[4] + e) % 2 ** 32
    '''for i in range(5):
        print(hex(h[i]))'''
    return h


def f(b, c, d, i):
    result = 0
    if 0 <= i < 20:
        result = (b & c) | ((~b) & d)
    elif 20 <= i < 40:
        result = b ^ c ^ d
    elif 40 <= i < 60:
        result = (b & c) | (b & d) | (c & d)
    elif 60 <= i < 80:
        result = b ^ c ^ d
    return result


def K(i):
    result = 0
    if 0 <= i < 20:
        result = 0x5a827999
    elif 20 <= i < 40:
        result = 0x6ed9eba1
    elif 40 <= i < 60:
        result = 0x8f1bbcdc
    elif 60 <= i < 80:
        result = 0xca62c1d6
    return result
def str_2_bin(strr):
    result = ''
    for i in strr:
        x = bytes(i, encoding="utf-8")
        for j in x:
            y = bin(j)[2:]
            while len(y) % 8 != 0:
                y = '0' + y
            result += y
    return result
def sign(q,a,M,x,k):
    t = str_2_bin(M)
    byte_num = len(t) // 8
    if len(t) % 8 != 0:
        byte_num += 1
    m = int("0x" + hashlib.sha256(int.to_bytes(int(t, 2), byte_num, byteorder='big')).hexdigest(), 16)
    s1 = gmpy2.powmod(a, k, q)
    s2 = gmpy2.invert(k, q - 1) * (m - x * s1) % (q - 1)
    return (s1, s2)
def Sign(q,a,M,x,k):
    t = str_2_bin(M)
    byte_num = len(t) // 8
    if len(t) % 8 != 0:
        byte_num += 1
    m = int("0x" + hashlib.sha256(int.to_bytes(int(t, 2), byte_num, byteorder='big')).hexdigest(), 16)
    s1 = gmpy2.powmod(a, k, q)
    s2 = gmpy2.invert(k, q - 1) * (m - x * s1) % (q - 1)
    return (s1, s2)
def Vrfy(q,a,M,y,s1,s2):
    t = str_2_bin(M)
    byte_num = len(t) // 8
    if len(t) % 8 != 0:
        byte_num += 1
    m = int("0x" + hashlib.sha256(int.to_bytes(int(t, 2), byte_num, byteorder='big')).hexdigest(), 16)
    v1 = gmpy2.powmod(a, m, q)
    v2 = gmpy2.powmod(y, s1, q) * gmpy2.powmod(s1, s2, q) % q
    if v1 == v2:
        return "True"
    else:
        return "False"
class sm4:
    def __init__(self, k, iv,mode):
        if len(hex(k))!=34:
            print("k的长度错误！")
        if len(hex(iv))!=34:
            print("初始向量iv的长度错误！")
        self.k = k
        self.iv = iv
        self.mode= mode
    def encrypt(self,M):
        if self.mode=="CBC":
            return CBC(self.k,M,self.iv,1)
        elif self.mode=="CTR":
            return CTR(self.k,M,self.iv,1)
        else:
            print("Error!")
    def decrypt(self,M):
        if self.mode=="CBC":
            return CBC(self.k,M,self.iv,0)
        elif self.mode=="CTR":
            return CTR(self.k,M,self.iv,0)
        else:
            print("Error!")
    def file_encrypt(self,inputpath,outputpath):
        if self.mode=="CBC":
            CBC_file_encrypt(inputpath,outputpath,self.k,self.iv)
        elif self.mode=="CTR":
            CTR_file_encrypt(inputpath,outputpath,self.k,self.iv)
        else:
            print("Error!")
    def file_decrypt(self,inputpath,outputpath):
        if self.mode=="CBC":
            CBC_file_decrypt(inputpath,outputpath,self.k,self.iv)
        elif self.mode=="CTR":
            CTR_file_decrypt(inputpath,outputpath,self.k,self.iv)
class RSA:
    def __init__(self, n):
        if n<=0:
            print("输入错误！")
        self.n=n
    def export_PrivateKey(self):
        self.p = getPrime(self.n)
        self.q = getPrime(self.n)
        return (self.p,self.q)
    def export_PublicKey(self):
        self.e = getPrime(self.n)
        (x, y, z) = egcd((self.p - 1) * (self.q - 1), self.e)
        self.d = int(y)
        if (self.d < 0):
            self.d = self.d + (self.p - 1) * (self.q - 1)
        return (self.e,self.q*self.p)
    def PrivateKey_import(self,p,q):
        if check(p)==False or check(q)==False or p<0 or q<0:
            print("输入错误：输入私钥不是素数")
        self.p = p
        self.q = q
    def PublicKey_import(self,e,N):
        if e<=0 or N<=0:
            print("输入错误！")
        self.e = e
        self.N = N
    def encrypt(self,s):
        plaintext = list(s.encode())
        plaintext_int = int("".join(["{:08b}".format(x) for x in plaintext]), 2)
        ciphertext = RSA_encrypt(plaintext_int,self.e,self.N)
        return ciphertext
    def decrypt(self,s):
        cliphertext_bin = bin(RSA_decrypt(s,self.e,self.p,self.q))[2:]
        if len(cliphertext_bin) % 8 != 0:
            cliphertext_bin = '0' * (8 - len(cliphertext_bin) % 8) + cliphertext_bin  # 补0
        plaintext = bytes(
            [int(cliphertext_bin[i * 8:(i + 1) * 8], 2) for i in range(len(cliphertext_bin) // 8)]).decode()
        return plaintext
class SHA1:
    def hash(self,m):
        return Hash(m)
    def hash_file(self,inputfile_path):
        bytes_list = bytes_read(inputfile_path)
        data = b""
        for q in bytes_list:
            data = data + q
        return(Hash1(hex(str2long(data))[2:]))
class ElGamal:
    def export_PrivateKey(self):
        n=random.randint(512,1024)
        self.p = getPrime(n)
        self.q = getPrime(160)
        while (self.p-1)%self.q!=0:
            print(self.q)
            self.q = getPrime(160)
        h=random.randint(1,self.p-1)
        self.g=gmpy2.powmod(h,(self.p-1)//self.q,self.p)
        return (self.p, self.q,self.g)
    def PrivateKey_import(self,q,a,x):
        if check(q)==False:
            print("输入错误：输入私钥不是素数")
        self.q=q
        self.a=a
        self.x=x
        self.y=pow(self.a,x,q)
    def PublicKey_import(self,q,a,y):
        self.q = q
        self.a = a
        self.x = y
    def sign(self,m):
        k=getPrime(1024)
        return Sign(self.q,self.a,m,self.x,k)
    def vrfy(self,m,s1,s2):
        return Vrfy(self.q,self.a,m,self.y,s1,s2)
if __name__ == "__main__":
    s = b"a80528bee168857bb6d2d5efaf9db759d2500bce1ab218391bdb0b839d5b2864e9531842fd0c386820b4d7def9c16868848ea11bab92465cde82ed1944947d32733c00b4919733d2940b9d6f0fb7bb05824a383d9ce77211e3ee838fb55a60c1bf162b151bd6ea202b18e06989b293a31319d90ead52fdc54e920b04b3cf50ed01b657d66d43874289ac4cca65ccd7d27194"
    print(CBC(0x7dfee2f716c4c4cd5217f0d57c75c2d7,s,0xa8638d2fb23cc49206edd7c84532eaab,0))