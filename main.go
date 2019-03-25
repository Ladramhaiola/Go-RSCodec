package main

import (
	"fmt"

	"./reedsolomon"
)

func main() {
	rs := &reedsolomon.RSCodec{
		Prim:       0x11d,
		FieldSize:  255,
		EccSymbols: 6,
	}
	rs.InitLookupTables()

	encoded := rs.Encode("hello Reed-Solomon world")
	fmt.Println("Encoded:        ", encoded)

	// corrupt the message
	encoded[0] = 20
	encoded[1] = 0
	encoded[2] = 3
	fmt.Println("Corrupted:      ", encoded)

	// decode & correct
	decoded, _ := rs.Decode(encoded)
	fmt.Println("Decoded bytes:  ", decoded)
	fmt.Print("Decoded message: ")
	for _, i := range decoded {
		fmt.Print(string(i))
	}
	fmt.Println()
}
