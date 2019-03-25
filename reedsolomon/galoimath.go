package reedsolomon

import "errors"

// ============================ //
//    Math in Galoi Fields      //
// ============================ //

func gfAddition(x, y int) int {
	return x ^ y
}

func gfSubstraction(x, y int) int {
	return x ^ y
}

func gfMultiplication(x, y int) int {
	if x == 0 || y == 0 {
		return 0
	}
	return exponents[logs[x]+logs[y]]
}

func gfDivision(x, y int) (int, error) {
	if y == 0 {
		return -1, errors.New("Zero division")
	}
	if x == 0 {
		return 0, nil
	}
	return exponents[(logs[x]+255-logs[y])%255], nil
}

func gfPow(x, power int) int {
	return exponents[negmod(logs[x]*power, 255)]
}

func gfInverse(x int) int {
	return exponents[255-logs[x]]
}

// multiply polynomial by scalar
func gfPolyScale(p []int, x int) []int {
	result := make([]int, len(p))
	for i := 0; i < len(p); i++ {
		result[i] = gfMultiplication(p[i], x)
	}
	return result
}

func gfPolyAddition(p, q []int) (result []int) {
	if len(p) > len(q) {
		result = make([]int, len(p))
	} else {
		result = make([]int, len(q))
	}
	for i := 0; i < len(p); i++ {
		result[i+len(result)-len(p)] = p[i]
	}
	for i := 0; i < len(q); i++ {
		result[i+len(result)-len(q)] ^= q[i]
	}
	return
}

// multiply two polynomials inside Galois Field
func gfPolyMultiplication(p, q []int) (result []int) {
	result = make([]int, len(p)+len(q)-1)
	for j := 0; j < len(q); j++ {
		for i := 0; i < len(p); i++ {
			result[i+j] ^= gfMultiplication(p[i], q[j])
		}
	}
	return
}

func gfPolyEvaluate(p []int, x int) int {
	// Evaluates a polynomial in GF(2^p) given the value for x.
	// This is based on Horner's scheme for maximum efficiency.
	// example: 01 x4 + 0f x3 + 36 x2 + 78 x + 40 = (((01 x + 0f) x + 36) x + 78) x + 40
	y := p[0]
	for i := 1; i < len(p); i++ {
		y = gfMultiplication(y, x) ^ p[i]
	}
	return y
}

func gfPolyDivision(divident, divisor []int) ([]int, []int) {
	// Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations
	result := make([]int, len(divident))
	copy(result, divident)

	for i := 0; i < len(divident)-(len(divisor)-1); i++ {
		coef := result[i]
		if coef != 0 {
			for j := 1; j < len(divisor); j++ {
				if divisor[j] != 0 {
					result[i+j] ^= gfMultiplication(divisor[j], coef)
				}
			}
		}
	}
	separator := len(divisor) - 1
	// return quotient, remainder
	return result[:len(result)-separator], result[len(result)-separator:]
}
