#!/usr/bin/python
import math
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

# find prime factors
def fast_primes(n):
	""" Returns  a list of primes < n """
	sieve = [True] * (n+1)
	for i in range(3,int(n**0.5)+1,2):
		if sieve[i]:
			#for multiples of 3, starts with 9, increase every 2x3
			sieve[i*i::2*i] = [False]*len(sieve[i*i::2*i])
	return [2]+[i for i in range(3,n+1,2) if sieve[i]]

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
print ('the first triangle number to have over five hundred divisors is ',Euler12(500))

'''   Question 13  '''
'''   Question 14  '''
'''   Question 15  '''
'''   Question 16  '''
'''   Question 17  '''
'''   Question 18  '''
'''   Question 19  '''
'''   Question 20  '''
'''   Question 21  '''
'''   Question 22  '''
'''   Question 23  '''
'''   Question 24  '''

#a = [("'''   Question %d  '''")%(i) for i in range(6,25,1)]
#for i in a:
#	print (i)