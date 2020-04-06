#!/usr/bin/python3
def get_cone_volume(radius, height):
	volume = (height/3)*3.1416*(radius**2)
	return volume

print(get_cone_volume(2,3))

def recursive_factorial(n):
	if n==1:
		return n
	else:
		return n*recursive_factorial(n-1)

print(recursive_factorial(5))


def factorial(m):
	num=1
	while m>1:
		number=m
		fac= num*number
		num=fac
		m=m-1
	return fac

print(factorial(5))

import time
def count_down(n, odd=False):
	if odd==False:
		while n>=0:
			print(n)
			n-=1
			time.sleep(0.1)
	else:
		while n>=0:
			if n%2==0:
				n-=1
			else:
				print(n)
				n-=1
			time.sleep(0.1)
count_down(20,True)

def get_final_price(price, discount_percentage=10):
	"Return the final price after applying the discount percentage"
	return (price  * (100-discount_percentage) / 100)
print(get_final_price(50))
