// r8brain-free-src Copyright (c) 2013-2014 Aleksey Vaneev
// See the "License.txt" file for license.
//
// Please read the "r8bsrc.h" file for function descriptions.

unit r8bsrc;

interface

type
	CR8BResampler = Pointer;
	PR8BDouble = ^Double;

const
	r8brr16 = 0; // 16-bit precision resampler.
	r8brr16IR = 1; // 16-bit precision resampler for impulse responses.
	r8brr24 = 2; // 24-bit precision resampler (including 32-bit floating
		// point).

function r8b_create( SrcSampleRate: Double; DstSampleRate: Double;
	MaxInLen: LongInt; ReqTransBand: Double; Res: LongInt ): CR8BResampler;
	cdecl; external 'r8bsrc.dll';

procedure r8b_delete( rs: CR8BResampler ); cdecl; external 'r8bsrc.dll';

function r8b_get_inlen( rs: CR8BResampler ): LongInt; cdecl;
	external 'r8bsrc.dll';

procedure r8b_clear( rs: CR8BResampler ); cdecl; external 'r8bsrc.dll';

function r8b_process( rs: CR8BResampler; ip0: PR8BDouble; l: LongInt;
	var op0: PR8BDouble ): LongInt; cdecl; external 'r8bsrc.dll';

implementation

end.
