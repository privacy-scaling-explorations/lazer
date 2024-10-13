package lazer

/*
#cgo CFLAGS: -g -Wall -Wextra -I../../
#cgo LDFLAGS: ./../liblazer.a -lmpfr -lgmp -lm

#include <gmp.h>
#include <mpfr.h>
#include "lazer.h"
*/
import "C"
import (
	"fmt"
	"unsafe"
)

const pubkeylen = 897
const privkeylen = 1281
const msglen = 64
const masked_msglen_max = 22000
const blindsiglen_max = 3000
const siglen_max = 32000

type UserStateT struct {
	c *C.user_state_struct
}

type SignerStateT struct {
	c *C.signer_state_struct
}

type VerifierStateT struct {
	c *C.verifier_state_struct
}

func UserInit(pubkey []byte) (state UserStateT) {

	_pubkey := C.CBytes(pubkey)

	if len(pubkey) != pubkeylen {
		panic(fmt.Sprintf("UserInit: len(pubkey) != %d", pubkeylen))
	}

	state.c = (*C.user_state_struct)(C.malloc(C.sizeof_user_state_struct))
	C.user_init(state.c, (*C.uint8_t)(_pubkey))

	C.free(_pubkey)
	return state
}

func UserMaskmsg(state *UserStateT, msg []byte) (masked_msg []byte) {
	_masked_msg := C.malloc(masked_msglen_max)
	_masked_msglen := C.malloc(C.sizeof_size_t)
	_msg := C.CBytes(msg)

	C.user_maskmsg(state.c, (*C.uint8_t)(_masked_msg), (*C.size_t)(_masked_msglen), (*C.uint8_t)(_msg))
	masked_msg = C.GoBytes(_masked_msg, (C.int)(*(*C.size_t)(_masked_msglen)))

	C.free(_masked_msg)
	C.free(_masked_msglen)
	C.free(_msg)
	return masked_msg
}

func UserSign(state *UserStateT, blindsig []byte) (rc int, sig []byte) {

	_sig := C.malloc(siglen_max)
	_siglen := C.malloc(C.sizeof_size_t)
	_blindsig := C.CBytes(blindsig)

	_rc := C.user_sign(state.c, (*C.uint8_t)(_sig), (*C.size_t)(_siglen), (*C.uint8_t)(_blindsig), 0)
	rc = int(_rc)
	if rc == 1 {
		sig = C.GoBytes(_sig, (C.int)(*(*C.size_t)(_siglen)))
	} else {
		sig = nil
	}

	C.free(_sig)
	C.free(_siglen)
	C.free(_blindsig)
	return rc, sig
}

func UserClear(state *UserStateT) {
	C.user_clear((*C.user_state_struct)(state.c))
	C.free((unsafe.Pointer)(state.c))
}

func SignerKeygen() (privkey []byte, pubkey []byte) {
	_pubkey := C.malloc(pubkeylen)
	_privkey := C.malloc(privkeylen)

	C.signer_keygen((*C.uint8_t)(_privkey), (*C.uint8_t)(_pubkey))
	privkey = C.GoBytes(_privkey, privkeylen)
	pubkey = C.GoBytes(_pubkey, pubkeylen)

	C.free(_privkey)
	C.free(_pubkey)
	return privkey, pubkey
}

func SignerInit(pubkey []byte, privkey []byte) (state SignerStateT) {

	_pubkey := C.CBytes(pubkey)
	_privkey := C.CBytes(privkey)

	if len(pubkey) != pubkeylen {
		panic(fmt.Sprintf("SignerInit: len(pubkey) != %d", pubkeylen))
	}
	if len(privkey) != privkeylen {
		panic(fmt.Sprintf("SignerInit: len(privkey) != %d", privkeylen))
	}

	state.c = (*C.signer_state_struct)(C.malloc(C.sizeof_signer_state_struct))
	C.signer_init(state.c, (*C.uint8_t)(_pubkey), (*C.uint8_t)(_privkey))

	C.free(_pubkey)
	C.free(_privkey)
	return state
}

func SignerSign(state *SignerStateT, masked_msg []byte) (rc int, blindsig []byte) {

	_blindsig := C.malloc(blindsiglen_max)
	_blindsiglen := C.malloc(C.sizeof_size_t)
	_masked_msg := C.CBytes(masked_msg)

	_rc := C.signer_sign(state.c, (*C.uint8_t)(_blindsig), (*C.size_t)(_blindsiglen), (*C.uint8_t)(_masked_msg), 0)
	rc = int(_rc)
	if rc == 1 {
		blindsig = C.GoBytes(_blindsig, (C.int)(*(*C.size_t)(_blindsiglen)))
	} else {
		blindsig = nil
	}

	C.free(_blindsig)
	C.free(_blindsiglen)
	C.free(_masked_msg)
	return rc, blindsig
}

func SignerClear(state *SignerStateT) {
	C.signer_clear(state.c)
	C.free((unsafe.Pointer)(state.c))
}

func VerifierInit(pubkey []byte) (state VerifierStateT) {

	_pubkey := C.CBytes(pubkey)

	if len(pubkey) != pubkeylen {
		panic(fmt.Sprintf("VerifierInit: len(pubkey) != %d", pubkeylen))
	}

	state.c = (*C.verifier_state_struct)(C.malloc(C.sizeof_verifier_state_struct))
	C.verifier_init(state.c, (*C.uint8_t)(_pubkey))

	C.free(_pubkey)
	return state
}

func VerifierVrfy(state *VerifierStateT, msg []byte, sig []byte) int {
	_msg := C.CBytes(msg)
	_sig := C.CBytes(sig)

	_rc := C.verifier_vrfy(state.c, (*C.uint8_t)(_msg), (*C.uint8_t)(_sig), 0)

	C.free(_msg)
	C.free(_sig)
	return int(_rc)
}

func VerifierClear(state *VerifierStateT) {
	C.verifier_clear(state.c)
	C.free((unsafe.Pointer)(state.c))
}
