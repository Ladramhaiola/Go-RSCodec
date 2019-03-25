package reedsolomon

import (
	"fmt"
)

// RSCodec Reed-Solomon coder/decoder
type RSCodec struct {
	// Primitive polynomial for lookup table generation
	Prim int
	// Fulle field size (eg 255 for 8-bit symbols)
	FieldSize int
	// Number of ECC symbols in message
	EccSymbols int
}

var (
	exponents = make([]int, 512)
	logs      = make([]int, 256)
)

// InitLookupTables fills exponential & log tables
func (r *RSCodec) InitLookupTables() {
	x := 1
	for i := 0; i < 255; i++ {
		exponents[i] = x
		logs[x] = i
		x = russianPeasantMult(x, 2, r.Prim, r.FieldSize+1, true)
	}

	for i := 255; i < 512; i++ {
		exponents[i] = exponents[i-255]
	}
}

// Encode given message into Reed-Solomon
func (r *RSCodec) Encode(data string) (encoded []int) {
	byteMessage := make([]int, len(data))
	for i, ch := range data {
		byteMessage[i] = int(ch)
	}
	fmt.Println("Original byte message:", byteMessage)

	generator := rsGeneratorPoly(r.EccSymbols)

	placeholder := make([]int, len(generator)-1)
	// Pad the message and divide it by the irreducible gnerator polynomial
	_, remainder := gfPolyDivision(append(byteMessage, placeholder...), generator)

	encoded = append(byteMessage, remainder...)
	return
}

// Decode and correct encoded Reed-Solomon message
func (r *RSCodec) Decode(data []int) ([]int, []int) {
	decoded := data

	synd := calcSyndromes(data, r.EccSymbols)
	if checkSyndromes(synd) {
		m := len(decoded) - r.EccSymbols
		return decoded[:m], decoded[m:]
	}

	fsynd := calcForneySyndromes(synd, []int{}, len(decoded))

	// compute the error locator polynomial using Berlekamp-Massey
	errLoc := unknownErrorLocator(fsynd, r.EccSymbols)
	// reverse errLoc
	for i, j := 0, len(errLoc)-1; i < j; i, j = i+1, j-1 {
		errLoc[i], errLoc[j] = errLoc[j], errLoc[i]
	}
	errPos := findErrors(errLoc, len(decoded))

	decoded = correctErrors(decoded, synd, errPos)
	synd = calcSyndromes(decoded, r.EccSymbols)

	m := len(decoded) - r.EccSymbols
	return decoded[:m], decoded[m:]
}

func rsGeneratorPoly(nsym int) []int {
	// generate an irreducible polynomial (necessary to encode message in Reed-Solomon)
	g := []int{1}
	for i := 0; i < nsym; i++ {
		g = gfPolyMultiplication(g, []int{1, gfPow(2, i)})
	}
	return g
}
