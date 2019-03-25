package reedsolomon

func calcSyndromes(message []int, nsym int) []int {
	synd := make([]int, nsym)
	for i := 0; i < nsym; i++ {
		synd[i] = gfPolyEvaluate(message, gfPow(2, i))
	}
	synd = append([]int{0}, synd...)
	return synd
}

func calcForneySyndromes(synd, pos []int, messageLen int) []int {
	// Compute Forney syndromes, which computes a modified syndromes to compute only errors (erasures are trimmed out).
	erasePosReversed := make([]int, len(pos))
	for i, p := range pos {
		// need to convert the positions to coefficients degrees for the errata locator algo to work (eg: instead of [0, 1, 2] it will become [len(msg)-1, len(msg)-2, len(msg) -3])
		erasePosReversed[i] = messageLen - 1 - p
	}

	fsynd := make([]int, len(synd)-1)
	copy(fsynd, synd[1:])

	for i := 0; i < len(pos); i++ {
		x := gfPow(2, erasePosReversed[i])
		for j := 0; j < len(fsynd); j++ {
			fsynd[j] = gfMultiplication(fsynd[j], x) ^ fsynd[j+1]
		}
	}
	return fsynd
}

func checkSyndromes(synd []int) bool {
	for _, v := range synd {
		if v > 0 {
			return false
		}
	}
	return true
}
