package reedsolomon

import (
	"fmt"
	"log"
)

// RSCodec Reed-Solomon coder/decoder
type RSCodec struct {
	// Primitive polynomial for lookup table generation
	Primitive int
	// Number of ECC symbols in message
	EccSymbols int
}

var (
	exponents = make([]int, 512)
	logs      = make([]int, 256)
)

// InitLookupTables fills exponential & log tables
func (r *RSCodec) InitLookupTables() {
	// Precompute the logarithm and anti-log tables for faster computation, using the provided primitive polynomial.
	// The idea: b**(log_b(x), log_b(y)) == x * y, where b is the base or generator of the logarithm =>
	// we can use any b to precompute logarithm and anti-log tables to use for multiplying two numbers x and y.
	x := 1
	for i := 0; i < 255; i++ {
		exponents[i] = x
		logs[x] = i
		x = russianPeasantMult(x, 2, r.Primitive, 256, true)
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
	fmt.Println("Original message:", byteMessage)

	g := rsGeneratorPoly(r.EccSymbols)

	placeholder := make([]int, len(g)-1)
	// Pad the message and divide it by the irreducible gnerator polynomial
	_, remainder := gfPolyDivision(append(byteMessage, placeholder...), g)

	encoded = append(byteMessage, remainder...)
	return
}

// Decode and correct encoded Reed-Solomon message
func (r *RSCodec) Decode(data []int) ([]int, []int) {
	decoded := data
	if len(data) > 255 {
		log.Fatalf("Message is too long, max allowed size is %d\n", 255)
	}

	synd := calcSyndromes(data, r.EccSymbols)
	if checkSyndromes(synd) {
		m := len(decoded) - r.EccSymbols
		return decoded[:m], decoded[m:]
	}

	// compute the error locator polynomial using Berlekamp-Massey
	errLoc := unknownErrorLocator(synd, r.EccSymbols)
	// reverse errLoc
	reverse(errLoc)
	errPos := findErrors(errLoc, len(decoded))

	decoded = correctErrors(decoded, synd, errPos)

	synd = calcSyndromes(decoded, r.EccSymbols)
	if !checkSyndromes(synd) {
		log.Fatalf("Could not correct message\n")
	}

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
