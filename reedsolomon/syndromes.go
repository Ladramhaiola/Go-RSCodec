package reedsolomon

func calcSyndromes(message []int, nsym int) []int {
	synd := make([]int, nsym)
	for i := 0; i < nsym; i++ {
		synd[i] = gfPolyEvaluate(message, gfPow(2, i))
	}
	synd = append([]int{0}, synd...)
	return synd
}

func checkSyndromes(synd []int) bool {
	for _, v := range synd {
		if v > 0 {
			return false
		}
	}
	return true
}
