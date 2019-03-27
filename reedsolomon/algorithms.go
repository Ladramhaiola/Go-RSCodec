package reedsolomon

import (
	"log"
	"math"
)

// ========================================== //
//             Used Algorithms                //
// ========================================== //

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

// rotate the array
func reverse(a []int) {
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
}

func filter(a []int, f func(int) bool) []int {
	b := make([]int, 0)
	for _, v := range a {
		if f(v) {
			b = append(b, v)
		}
	}
	return b
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
	// Forney algorithm to compute the magnitudes
	// E will store the values that need to be corrected (error magnitude) to correct the input message
	E := make([]int, len(message))

	for i, location := range locationPolynomial {
		locationInverse := gfInverse(location)

		// Compute the formal derivative of the error locator polynomial
		// the formal derivative of the locator is used as the denominator for Forney algorithm,
		// which simply says that the ith error value is given by error_evaluator(gf_inverse(Xi)/error_locator_deriv(gf_inverse(Xi)))
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

		if errorLocatorPrime == 0 {
			log.Println("Could not find magnitude")
		}

		// compute the magnitude
		// magnitude is the correction vector
		magnitude, _ := gfDivision(y, errorLocatorPrime)
		E[errPos[i]] = magnitude
	}

	return E
}

// FindPrimePolys computes the list of prime polynomials for the given generator
// and galois field characteristic exponent.
func FindPrimePolys(cExponent int, fast, single bool) []int {
	// Prime irreducible polynomial is used to avoid duplicate values in LUT
	fieldCharac := int(math.Pow(2.0, float64(cExponent))) - 1
	fieldCharacNext := int(math.Pow(2.0, float64(cExponent+1))) - 1

	primCandidates := []int{}
	if fast {
		// generate maybe prime polynomials and check later if they really are irreducible
		primCandidates = sieveOfEratosthenes(fieldCharacNext)
		primCandidates = filter(primCandidates, func(x int) bool { return x > fieldCharac })
	} else {
		// try each possible prime polynomial
		for i := fieldCharac + 2; i < fieldCharacNext; i += 2 {
			primCandidates = append(primCandidates, i)
		}
	}

	correctPrimes := []int{}
	for _, prim := range primCandidates {
		seen := make([]int, fieldCharac+1)
		conflict := false

		// second loop, build the whiole Galoi field
		x := 1
		for i := 0; i < fieldCharac; i++ {
			// compute the next value in the field (ie, the next power of the generator)
			x = russianPeasantMult(x, 2, prim, fieldCharac+1, true)

			// Rejection criterion: if the value overflowed (above fieldCharac) or is duplicate of a
			// previously generated power of alpha, then we reject the polynomial (not prime)
			if x > fieldCharac || seen[x] == 1 {
				conflict = true
				break
			} else {
				seen[x] = 1
			}
		}

		if !conflict {
			correctPrimes = append(correctPrimes, prim)
			if single {
				return []int{prim}
			}
		}
	}

	return correctPrimes
}
