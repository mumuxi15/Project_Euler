#!/usr/bin/env python3

#!/usr/bin/python
import math

###  ---------------- HELPER FUNCTIONS -------------------   ###
def permutation(lst):
    if len(lst)==1:
        print ('return : ',lst)
        return [lst]
    if len(lst) == 0:
        return []
    
    
    """  [A|B C], [B|A,C]; fixed + recursive"""
    l = []
    for i in range(len(lst)):
        fixed = lst[i]
        r = lst[0:i]+lst[i+1:]
#		print (fixed,' | ',r)
        for p in permutation(r):
            l.append([fixed]+p)
            
#	print ('-----'*10,l)
    return l



# find prime factors
def fast_primes(n):
    """ Returns  a list of primes < n """
    sieve = [True] * (n+1)
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i]:
            #for multiples of 3, starts with 9, increase every 2x3
            sieve[i*i::2*i] = [False]*len(sieve[i*i::2*i])
    return [2]+[i for i in range(3,n+1,2) if sieve[i]]

###  ----------------------------------   ###



'''   Question 1  '''

#def Euler1(N):
#	"""  failed at larger N due to large rounding error in python3 """
#	m3 = [i*3 for i in range(int((N-0.1)//3)+1)]
#	m5 = [i*5 for i in range(int((N-0.1)//5)+1)]
#	rs = set(m3+m5)
#	return sum(rs)

def Euler1(N):
    " algorithmic approach - all pass"
    m3 = int((N-0.1)//3)
    m5 = int((N-0.1)//5)
    m15 = int((N-0.1)//15)
    m3_sum = 3*(1+m3)*m3//2
    m5_sum = 5*(1+m5)*m5//2
    m15_sum = 15*(1+m15)*m15//2
    rs = m3_sum+m5_sum-m15_sum
    return rs

#print ('1. Sum of all the multiples of 3 or 5 below 1000: ',Euler1(1000))


'''   Question 2  Fibonacci sequence  '''
def Euler2(x):
    """  F(1)=0, F(2)=1  
        F(n) = F(n-1)+F(n-2)  """
    a, b = 0, 1
    even = 0
    l = [a] # fib list
    while b <= x:
        a, b = b, a+b
        if a%2 ==0:
            even+=a
    return even

#print ('2.Sum of even term of fib series under 4 million: ',Euler2(400000000))

'''   Question 3 Largest Prime '''
def LargestPrime(n):
    ''' start with largest prime number
    even number remove 2 until odd
    '''
    maxPrime = 1
    while n % 2 == 0:
        n = n//2
        maxPrime = 2
    for i in range(3,int(n**0.5)+1,1):
        while n % i ==0: 
            n = n//i
            maxPrime = i
    if n > 2:
        maxPrime  = n
    return maxPrime
#print ('3.Largest prime factor of the number 600851475143: ',LargestPrime(600851475143))

'''   Question 4  '''

def panlidromic(x):
    maxP = 1
    for i in range(int(x**0.5),100,-1):
        for j in range(min(x//i,999),100,-1): #3 digits
            p = i*j
            if p > maxP:
                s = str(p)
                if s == s[::-1]:
                    maxP = max(maxP, p)
    return None if maxP <=1 else maxP
#  ----------------------------------
#  quick method:
#  palindromic = 1e5*a + 1e4*b + 1e3*c + 1e2*b + 10a 
#              = 100001a + 10010b + 1100 c 
#              = 11(9091a+910b+100c)
#  ----------------------------------

def panlidromic2():
    # check if i or j is divisible by 11
    maxP = 1
    for i in range(990,110,-11):
        for j in range(999,100,-1): #3 digits
            p = i*j
            if p > maxP:
                s = str(p)
                if s == s[::-1]:
                    maxP = max(maxP, p)
    return None if maxP <=1 else maxP

#print ('4. The largest palindrome made from the product of two 3-digit numbers: ',panlidromic2())

'''   Question 5  '''
'''2520 is the smallest number that can be divided by each of the numbers from 1 to 10 without any remainder.

What is the smallest positive number that is evenly divisible by all of the numbers from 1 to 20?
'''

def smallest_divider(x):
    '''find all primes <= x '''
    prime = {x:1 for x in fast_primes(x)}
    prime[2] = int(math.log(x,2))
    prime_list = [x for x in prime.keys() if x != 2]
    for u in range(3,x+1,2):
        for p in prime_list:
            print ('--',u,p)
            if u>p:
                print (u,p)
                c = 0
                while u % p == 0:
                    u = u//p
                    c += 1
                prime[p] = max(c,prime[p])
                print (prime)
                
    product = 1
    for k, v in prime.items():
        product = k**v*product
    return product

# gcd method 

def gcd(x,y):
    """ explanation: Every common divisor of a and b also div (b-a)
    """
    while y!=0:
        x,y=y,x%y
    return x

def smallest_divider(x):
    '''p = A X B / gcd(A,B)'''
    p = 2**int(math.log(x,2))
    for a in range(3,x+1,1):
        p = a*p//gcd(a,p)
    return p


#print ('5. The smallest positive number that is evenly divisible by all of the numbers from 1 to 20 ',smallest_divider(20))

'''   Question 6  '''
def Euler6(N):
    return (N*(N+1)//2)**2-N*(N+1)*(2*N+1)//6

'''   Question 7  '''
def Euler7(x):
    scale_factor = 3
    prime_list = []
    count = 0
    while len(prime_list)<x:
        prime_list = fast_primes(x*scale_factor)
        scale_factor = scale_factor *10
        count +=1
    return prime_list[x]
#print ('What is the 10001st prime number? ',Euler7(10001))

'''   Question 8  '''
def Euler8(x,k):
    s = [int(i) for i in str(x)]
    maxp = -1
    for i in range(len(s)-k+1):
        p = 1
        if 0 in s[i:i+k]:
            p = 0
        else:
            for j in s[i:i+k]:
                p = p*j
        print (s[i:i+k], ' product = ',p)
        maxp = max(p,maxp)
        
        
        '''   Question 9  '''
def pythagorean(N):
    '''a < N/3
        b > a 
        a+b > c > N/3
        b > N/6
    
    '''
    for a in range(N//3,0,-1):
        for b in range(a+1,(N-a)//2+1,1):
            c = N-a-b
            if (c > b)&( b+a>c):
                if a**2+b**2==c**2:
                    print ('triangle',a,b,c)
                    return a*b*c
                
    return -1

#print ('Find the product abc of triangle for which a+b+c = 70', pythagorean(70))

'''   Question 10  '''
def Euler10(n):
    """ Returns  a list of primes < n """
    sieve = [True] * (n+1)
    sieve[4::2] = [False] * (n//2-1)
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i]:
            #for multiples of 3, starts with 9, increase every 2x3
            sieve[i*i::2*i] = [False]*len(sieve[i*i::2*i])
    summ = 2
    l ={2:2}
    
    for x in range(3,n):
        if sieve[x]:
            summ += x 
        l[x] = summ
    return l

#print ('Find the sum of all the primes below two million: ',Euler10(2000000)[2000000])
'''   Question 11  '''
def Euler11(A):
    N = len(A)
    maxP = 0
    for j in range(N):
        for i in range(N):
            
            if i+4<=N:
                H = [A[j][i],A[j][i+1],A[j][i+2],A[j][i+3]]
                if 0 in H:
                    continue 
                
                maxP = max(maxP, A[j][i]*A[j][i+1]*A[j][i+2]*A[j][i+3])
                #					print (H, ' s= ',maxS, 'p = ',maxP//10000)
                
            if j+4<=N:
                V = [A[j][i],A[j+1][i],A[j+2][i],A[j+3][i]]
                if 0 in V:
                    continue 
                maxP = max(maxP, A[j][i]*A[j+1][i]*A[j+2][i]*A[j+3][i])
                
            if (i+4<=N)&(j+4<=N):
                D = [A[j][i],A[j+1][i+1],A[j+2][i+2],A[j+3][i+3]]
                if 0 in D:
                    continue 
                maxP = max(maxP, A[j][i]*A[j+1][i+1]*A[j+2][i+2]*A[j+3][i+3])
                
            if (j>=3)&(i+4<=N):
                D= [A[j][i],A[j-1][i+1],A[j-2][i+2],A[j-3][i+3]]
                if 0 in D:
                    continue 
                maxP = max(maxP,A[j][i]*A[j-1][i+1]*A[j-2][i+2]*A[j-3][i+3])
                
    return maxP

'''   Question 12  '''

def find_factors(x):
    factors = 1
    for i in range(2,int(math.sqrt(x))+1,1):
        if x % i ==0:
            factors += 2 if i != x//i else 1
    return factors

def Euler12(n):
    l = {}
    for x in range(2,50000,1):
        l[x] = find_factors(x)
    l[1] = 0
    
    factors = 0
    i = 2
    while factors <n:
        if i%2 ==0:
            factors = (l[i//2]+1)*(l[i+1]+1)-1
        else:
            factors = (l[i]+1)*(l[(i+1)//2]+1)-1
        i += 1
    return i*(i-1)//2
#print ('the first triangle number to have over five hundred divisors is ',Euler12(500))

'''   Question 13  '''
def Euler13():
    s = sum([int(input()) for l in range(int(input()))])
    print (str(s)[0:10])
    
    '''   Question 14  '''
    
def collatz(n):
    return n//2 if n%2==0 else 3*n+1


"""define list of steps to reach 1 for all numbers are 0
    for each number n, 
    record all the intermediates and store it in a list
    if x > N, do not skip recording until x falls under N
    calculate steps for each intermediates
    compare maxchain and maxid, update if new high
"""
def Euler14(N):
    l = [0]*(N+1)
    l[2] = 1
    maxid_list = [0,1]
    maxid, maxchain = 1,0
    
    for n in range(2,N+1,1):
        c, x = 0, n
        steps = [[c,x]]
        while l[x] == 0:
            x = collatz(x)
            c += 1
            while x>N:
                x = collatz(x)
                c += 1
            steps.append([c,x])
            
            
        if len(steps)>1:
            for s, v in steps[:-1]:
                l[v] = l[steps[-1][-1]]+steps[-1][0]-s
        if maxchain <= l[n]:  #x is the last number of steps with max counts
            maxchain, maxid = l[n], n
#           print ('n: %d |  maxchain:%d|maxid: %d'%(n,maxchain,maxid))
#           print (steps,'\n')
        maxid_list.append(maxid)
    return maxid_list
#print ('Under one million, %d produces the longest chain'%(Euler14(1000000)[-1]))


'''   Question 15  '''

def print_grid(grid):
    for _ in grid:
        print (_)
    print ('-'*20)
    
def euler15(N,M):
    #   f(i,j) = f(i-1,j)+f(i,j-1)
    grid = [[1]*(N+1)]+[[1]+[0]*N for i in range(M)]
    for j in range(1,M+1):
        for i in range(1,N+1):
            grid[j][i] = grid[j-1][i]+grid[j][i-1]
    return grid[M][N]

# -----------   quick method  -----------  #
"""  it takes 40 steps to reach to the destination, and out of them 20 are arrow right. So the question becomes choosing 20 out of 40. 20C40 """
def fac(x):
    rs = 1
    for i in range(x,1,-1):
        rs = i*rs
    return rs

def Euler15_quick(N,M):
    
    MOD = 10**9+7
    return fac(M+N)//fac(M)//fac(N)%MOD
#print ('How many such routes are there through a 20×20 grid? ',euler15(20,20))

'''   Question 16  '''
#print (sum([int(a) for a in str(2**1000)])
    
'''   Question 17  '''
class Euler17:
    dict_a ={0:'',1:'One',2:'Two',3:'Three',4:'Four',5:'Five',6:'Six',7:'Seven',8:'Eight',9:'Nine',10:'Ten',11:'Eleven',12:'Twelve',13:'Thirteen',14:'Fourteen',15:'Fifteen', 16:'Sixteen', 17:'Seventeen',
        18:'Eighteen',19:'Eighteen'}
    dict_b={2:'Twenty',3:'Thirty',4:'Forty',5:'Fifty',6:'Sixty',7:'Seventy',8:'Eighty',9:'Ninety'}
    
    def read_hundred(self,x):
        word = []
        c,b,a = (x%1000)//100,(x%100)//10,x%10
        if c>0:
            word += [self.dict_a[c],'Hundred']
        if b>=2:
            word += [self.dict_b[(x%100)//10], self.dict_a[x%10]]
        else:
            word += [self.dict_a[x%100]]
        return word
    
    def to_words(self,x):
        words = []
        if x >= 1000000000:
            words += self.read_hundred(x//1000000000) +['Billion']
            x = x%1000000000
        if x >= 1000000:
            words += self.read_hundred((x%1000000000)//1000000)+ ['Million']
            x = x%1000000
        if x >= 1000:
            words += self.read_hundred((x%1000000)//1000)+ ['Thousand']
            x = x%1000
        if x<1000:
            words += self.read_hundred(x) 
        words = ' '.join([w for w in words if w !=''])
        print (words)
        
#E17 = Euler17()
#E17.to_words(1100020098)

'''   Question 18  '''
def Euler18(A):
    n = len(A)
    for i in range(n-2,-1,-1):
        for j in range(len(A[i])):
            A[i][j] = max(A[i][j]+A[i+1][j],A[i][j]+A[i+1][j+1])
    print (A[0][0])
    
    
'''   Question 19  '''

def Euler19(y1, m1, d1, y2, m2, d2):
    #1 Jan 1900 is Monday
    y0 = 1900
    duration = 2800     #weekday patterns repeat every 2800 years
#   months=np.zeros(duration*12)
#   mark =np.zeros(duration*12)
#   months[0::2]=31    #Jan
#   months[1::12]=28   #Feb
#   months[3::12]=30   #Apr
#   months[5::12]=30   #Jun
#   months[6::12]=31   #Jul
#   months[7::12]=31   #Aug
#   months[8::12]=30   #Sep
#   months[9::12]=31   #Oct
#   months[10::12]=30
#   months[11::12]=31
#   array of days in months 
    months=[0]*duration*12
    mark =[0]*duration*12
    months[0::2] =[31]*len(months[0::2])    #Jan
    months[1::12]=[28]*len(months[1::12])   #Feb
    months[3::12]=[30]*len(months[3::12])   #Apr
    months[5::12]=[30]*len(months[5::12])   #Jun
    months[6::12]=[31]*len(months[6::12])   #Jul
    months[7::12]=[31]*len(months[7::12])   #Aug
    months[8::12]=[30]*len(months[8::12])   #Sep
    months[9::12]=[31]*len(months[9::12])   #Oct
    months[10::12]=[30]*len(months[10::12])
    months[11::12]=[31]*len(months[11::12])
    # 1992 is leap year. in x years,  
    for i in range(0,duration):
        # if leap year Feb is 29
        if (y0+i)%4 == 0 and ((y0+i)%100!=0 or (y0+i)%400==0):
            months[12*i+1] = 29
            
    sum_days = 0
    for i in range(len(months)):
        sum_days+=months[i]
        if sum_days%7==6:
            mark[i+1] = 1 #sum of Jan-March is Apr 1st
    mark = mark*5 # backup twice is enough, but just in case
    
    
    #the weekday calendar repeat every 2800 years. creates 4816 sundays
    multiple = lambda y: (y-1900)//2800
    
    multiple1 = multiple(y1)
    multiple2 = multiple(y2)
    y1 = y1 - multiple1*2800
    y2 = y2 - multiple2*2800
    
    while y2<y1:
        multiple2 += -1
        y2 += 2800
        
        
    sundays = (multiple2-multiple1)*4816
    
    ## calculaltion within the array length
    i,j = (y1-y0)*12+m1-1,(y2-y0)*12+m2
    if d1>1:  #01/03 = start on 02/01
        i+=1
        
#   print (y1,m1,d1,' index ',i, '   ', y2,m2,d2,' index ',j,' sum:', sum(mark[i:j]))
    if len(mark[i:j])>0:
        return (sundays+sum(mark[i:j]))
    else:
        return (sundays+0)


#print ('How many Sundays fell on the first of the month during the twentieth century (1 Jan 1901 to 31 Dec 2000)? ', Euler19(1901 ,1, 1, 2000, 12, 31) )

'''   Question 21  '''
def find_all_factors(x):
    factors=[]
    for i in range(2,int(math.sqrt(x))+1,1):
        if x % i ==0:
            factors += [i,x//i]
    return set(factors)

def Euler21(N):
    maximum = N
    records = [0]*N
    amicable = []
    for i in range(3,maximum):
        sum_factors = 1+sum(find_all_factors(i))
        
        if sum_factors>i:  #only record sum > its own
            records[i] = sum_factors
#           print (i, '  - - -   record - - ')
        elif sum_factors<i: #check records
            if records[sum_factors] ==i:
                amicable+=[i,sum_factors]
#               print ('sum <  i  ',records[sum_factors])
    return amicable

#print ('The sum of all the amicable numbers under 10000 is ', sum(Euler21(10000)))

'''   Question 22  '''
def mapping_func22(name):
    return sum([ord(x.lower())-96 for x in name])

def Euler22(names):
    names = sorted(names)
    score = 0
    for i in range(len(names)):
        score += mapping_func22(names[i])*(1+names.index(names[i]))
    print (score)
    
#Euler22(names = ["ZO","SHON","LYNWOOD","JERE","HAI","ELDEN","DORSEY","DARELL","BRODERICK","ALONSO"])
    
    
'''   Question 23  '''
    
def sum_of_divisors(x):
    s = 0
    for i in range(2,int(math.sqrt(x))+1,1):
        if x % i ==0:
            s += i if i==x//i else x//i+i
    return s+1

def Euler23(N=28123):
    # all numbers > 28123 can be written as the sum of two abundant numbers    
    abundant = [i for i in range(3,N) if sum_of_divisors(i)>i ]
        
    is_abundant_sum = [False]*(N)
    for i in abundant:
        for j in abundant:
            z = i+j
            if z <N and is_abundant_sum[z] is False:
                is_abundant_sum[z]=True
                
    return sum([i for i, x in enumerate(is_abundant_sum) if x is False])
    
#print ('The sum of all the positive integers which cannot be written as the sum of two abundant number is ', Euler23(N=28123))
    
'''   Question 24  '''
    
def lexicographic_permutations(s,n):
    N = len(s)
    total = math.factorial(N)
    print ('total permutations: ',total)
    l = []
    
    for i in range(N,1,-1):
        idx = n//math.factorial(i-1)
        rem = n%math.factorial(i-1)
        print ('i=',i,idx,rem)
        print (math.factorial(i-1))

#lexicographic_permutations(s="0123",n=5)


'''   Question 25  '''
def Euler25(n):
    sys.set_int_max_str_digits(5001)
    a = 1
    b = 1
    count = 2
    rs = [0,1]
    length = 1
    c = 10
    
    for count in range(2,25000):
        if b>c:
            rs.append(count)
            c *= 10
        a, b = b, a + b  # Update Fibonacci numbers
    
    return rs[n]

#print ('The index of the first term in the Fibonacci sequence to contain 1000 digits is ', Euler25(1000))

'''  Question 26   '''

def get_prime(n):
    # find the greatest prime under n
    sieve = [True]*n
    for i in range(3,len(sieve),2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*len(sieve[i*i::2*i])
    sieve[4::2] = [False]*len(sieve[4::2])
    
    return sieve

def Euler26(n):
    ''' the longest recur would be the close reptend prime - 
    A prime p for which 1/p has a maximal period decimal expansion of p-1 digits
    '''
    #pow(4, 3, 5)  = 4^3%5
    n = n-1 # d<N 
    if n<8:
        return 3
    for i in range(n,2,-1):
        if prime_list[i]:
            recur = 1
            k = 1
            while pow(10,k,i)!=1:
                recur+=1
                k+=1
            if recur+1 == i:
                return i
            
#prime_list = get_prime(10001)
#print ('Find the value of d<1000 for which 1/d contains the longest recurring cycle in its decimal fraction part, d=', Euler26(1000))


'''  Question 27   '''
            