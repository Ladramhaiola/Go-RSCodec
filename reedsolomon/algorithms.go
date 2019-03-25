package reedsolomon

// efficient implementation of Sieve algo
// return list of primes less then N
func sieveOfEratosthenes(N int) (primes []int) {
	b := make([]bool, N)
	for i := 2; i < N; i++ {
		if b[i] == true {
			continue
		}
		primes = append(primes, i)
		for k := i * i; k < N; k += i {
			b[k] = true
		}
	}
	return
}

// python-like modular division
func negmod(d, m int) int {
	var res = d % m
	if (res < 0 && m > 0) || (res > 0 && m < 0) {
		return res + m
	}
	return res
}

// Russian Peasant Multiplicaton algorithm in GF
// default values for 8-bit size: fieldsize=256, carryless=true
func russianPeasantMult(x, y, prim, fieldsize int, carryless bool) (result int) {
	for y > 0 {
		if (y & 1) != 0 {
			// y is odd, then the corresponding x to y
			// (the sum of all x's corresponding to odd y's will give the final product)
			if carryless {
				result ^= x
			} else {
				result += x
			}
		}
		y >>= 1 // y mod 2
		x <<= 1 // x power 2
		// GF modulo: if x >= 256 then apply modular reduction using the primitive polynomial
		// if primitive number > 256 directly XOR
		if prim > 0 && x&fieldsize != 0 {
			x ^= prim
		}
	}
	return
}

func forney(message, errorPolynomial, locationPolynomial, errPos []int) []int {
	E := make([]int, len(message))

	for i, location := range locationPolynomial {
		locationInverse := gfInverse(location)

		// Compute the formal derivative of the error locator polynomial
		errorLocatorPrimeTemp := []int{}

		for j := 0; j < len(locationPolynomial); j++ {
			if j != i {
				errorLocatorPrimeTemp = append(errorLocatorPrimeTemp, gfSubstraction(1, gfMultiplication(locationInverse, locationPolynomial[j])))
			}
		}

		// compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
		errorLocatorPrime := 1

		for _, coef := range errorLocatorPrimeTemp {
			errorLocatorPrime = gfMultiplication(errorLocatorPrime, coef)
		}

		// Y1 = omega(X1.inverse()) / prod(1 - Xj*X1.inverse()) for j in len(X)
		y := gfPolyEvaluate(errorPolynomial, locationInverse)
		y = gfMultiplication(gfPow(location, 1), y)

		// compute the magnitude
		magnitude, _ := gfDivision(y, errorLocatorPrime)
		E[errPos[i]] = magnitude
	}

	return E
}
