package reedsolomon

import "log"

// compute the errors locator polynomial from the errors positions as input
func calcErrorLocatorPoly(errorPositions []int) []int {
	erasureLocations := []int{1}
	// erasures location = product(1 - x*alpha**i) for i in error positions (alpha is the alpha choosen to eval polynomials)
	for _, p := range errorPositions {
		erasureLocations = gfPolyMultiplication(erasureLocations, gfPolyAddition([]int{1}, []int{gfPow(2, p), 0}))
	}
	return erasureLocations
}

// compute the error evaluator polynomial Omega from the
// syndrome locator Sigma
func calcErrorPoly(synd, erasureLocations []int, nsym int) []int {
	// Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
	placeholder := make([]int, nsym+1)
	placeholder = append([]int{1}, placeholder...)

	_, remainder := gfPolyDivision(gfPolyMultiplication(synd, erasureLocations), placeholder)

	return remainder
}

// Find error locator and evaluator polynomials with Berlekamp-Massey algorithm
func unknownErrorLocator(synd []int, nsym int) []int {
	// The idea is that BM will iteratively estimate the error locator polynomial
	// To do this, it will compute a Discrepancy term called Delta, which will tell
	// if the error locator polynomial needs update or not
	errLoc := []int{1} // Sigma
	oldLoc := []int{1} // BM is an iterative algorithm, prev iteration values of Sigma

	syndShift := 0
	if len(synd) > nsym {
		syndShift = len(synd) - nsym
	}

	for i := 0; i < len(synd); i++ {
		K := i + syndShift
		delta := synd[K]
		for j := 1; j < len(errLoc); j++ {
			delta ^= gfMultiplication(errLoc[len(errLoc)-(j+1)], synd[K-j])
		}

		// Shift polynomials to compute next degree
		oldLoc = append(oldLoc, 0)

		// iteratively estimate the errata locator and evaluator polynomials
		if delta != 0 {
			if len(oldLoc) > len(errLoc) {
				// computing Sigma
				newLoc := gfPolyScale(oldLoc, delta)
				oldLoc = gfPolyScale(errLoc, gfInverse(delta))
				errLoc = newLoc
			}

			// update with the discrepancy
			errLoc = gfPolyAddition(errLoc, gfPolyScale(oldLoc, delta))
		}
	}

	// drop leading zeroes
	for len(errLoc) > 0 && errLoc[0] == 0 {
		errLoc = errLoc[1:]
	}

	errs := len(errLoc) - 1
	if (errs * 2) > nsym {
		log.Printf("Too many errors to correct: %d\n", errs)
	}

	return errLoc
}

func findErrors(errLoc []int, messageLen int) []int {
	errs := len(errLoc) - 1
	errPos := []int{}

	for i := 0; i < messageLen; i++ {
		if gfPolyEvaluate(errLoc, gfPow(2, i)) == 0 {
			errPos = append(errPos, messageLen-1-i)
		}
	}

	if len(errPos) != errs {
		log.Println("Too many (or few) errors found by Chien Search")
	}
	return errPos
}

func correctErrors(message, synd, errPos []int) []int {
	coefPos := make([]int, len(errPos))

	for i, p := range errPos {
		coefPos[i] = len(message) - 1 - p
	}

	errorLocatorPolynomial := calcErrorLocatorPoly(coefPos)

	// reverse errLoc
	for i, j := 0, len(synd)-1; i < j; i, j = i+1, j-1 {
		synd[i], synd[j] = synd[j], synd[i]
	}
	errorPolynomial := calcErrorPoly(synd, errorLocatorPolynomial, len(errorLocatorPolynomial)-1)

	locationPolynomial := []int{}
	for i := 0; i < len(coefPos); i++ {
		l := 255 - coefPos[i]
		locationPolynomial = append(locationPolynomial, gfPow(2, -l))
	}

	// Forney algorithm: compute the magnitudes
	E := forney(message, errorPolynomial, locationPolynomial, errPos)

	// Simply add correction vector to our message
	message = gfPolyAddition(message, E)
	return message
}
