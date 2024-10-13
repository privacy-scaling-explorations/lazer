package main

import (
	"crypto/rand"
	"fmt"
	"os"

	"XXX1.org/lazer"
)

func main() {
	iterations := 1

	fmt.Print("lazer blind-signature demo\n")
	fmt.Print("--------------------------\n\n")

	/* generate random message for demo */
	message := make([]byte, 64)
	_, err := rand.Read(message)
	if err != nil {
		panic("rand.Read failed")
	}

	fmt.Print("Signer generates keypair ... \n")
	privkey, pubkey := lazer.SignerKeygen()
	fmt.Print("[OK]\n\n")

	fmt.Print("Initialize user with public key ... \n")
	user := lazer.UserInit(pubkey)
	fmt.Print("[OK]\n\n")

	fmt.Print("Initialize signer with public and private key ... ")
	signer := lazer.SignerInit(pubkey, privkey)
	fmt.Print("[OK]\n\n")

	fmt.Print("Initialize verifier with public key ... ")
	verifier := lazer.VerifierInit(pubkey)
	fmt.Print("[OK]\n\n")

	for i := 0; i < iterations; i++ {
		fmt.Print("User outputs masked message (including a proof of its well-formedness) ... ")
		masked_msg := lazer.UserMaskmsg(&user, message)
		fmt.Print("[OK]\n")

		fmt.Printf("masked message (t,P1): %d bytes\n\n", len(masked_msg))

		fmt.Print("Signer checks the proof and if it verifies outputs a blind signature ... ")
		rc, blindsig := lazer.SignerSign(&signer, masked_msg)
		if rc != 1 {
			fmt.Print("masked message is invalid.\n")
			os.Exit(1)
		}
		fmt.Print("[OK]\n")

		fmt.Printf("blind signature (tau,s1,s2): %d bytes\n\n", len(blindsig))

		fmt.Print("User outputs a signature on the message ... ")
		rc, sig := lazer.UserSign(&user, blindsig)
		if rc != 1 {
			fmt.Print("decoding blindsig failed.\n")
			os.Exit(1)
		}
		fmt.Print("[OK]\n")

		fmt.Printf("signature (P2): %d bytes\n\n", len(sig))

		fmt.Print("Verfifier verifies the signature on the message ... ")
		rc = lazer.VerifierVrfy(&verifier, message, sig)
		if rc != 1 {
			fmt.Print("verification failed.\n")
			os.Exit(1)
		}
		fmt.Print("[OK]\n\n")
	}

	lazer.UserClear(&user)
	lazer.SignerClear(&signer)
	lazer.VerifierClear(&verifier)
	os.Exit(0)
}
