
#################################################################################
# Hashing based on operations (>>) and (<<) in the symmetric group
#################################################################################  

from hash_aux import ut, en, op

class hash_sn:

    def __init__(self, string, blocksize, factorial, prime):
        self.m = string
        self.s = blocksize // 8
        self.t = factorial
        self.p = prime
        self.len_orig = len(string)
        self.pad = None
        self.h = None
    
    # Padding
    def hash_sn_pad(self):

        s = self.s
        m = self.m
        len_orig = self.len_orig
        s_2 = s//2

        if len_orig % s < s_2: len_pad = s_2 - (len_orig % s)
        else: len_pad = s - (len_orig % s_2)

        if len_orig > s_2: len_min = s_2
        else:
            len_min = len_orig
            len_pad += s_2 - len_orig

        m += '\x80'
        for i in range(1, len_pad): m += '\x00'

        r=''
        for i in range(len_min):
            m_char = m[len_orig-1-i]
            m += m_char
            r += format(ord(m_char),'08b')

        # Reverse last half block rotated s/2 positions to the left
        r= en.bits_to_int(r) << (s_2 * 8)

        self.pad = m
        self.r = r

    def hash(self):

        hash_sn.hash_sn_pad(self)

        m = self.pad
        len_orig = self.len_orig
        t = self.t
        s = self.s
        r = self.r
        p = self.p
        s_2 = s // 2
        two_exp_s = 2**(s * 8)
        l = len_orig << (s_2 * 8)
        len_new = len(m)

        # Initial state: h = r << l
        h = op.a_rotate_b_inv(l, r, t)
        k = 1

        for i in range(len_new // s):
            ind = i*s
            B = en.string_to_bits(m[ind : ind + s])
            B = en.bits_to_int(B)

            # Accumulative product of blocks modulo p
            if B: k = (k * B) % p
            else:  k = (k * k) % p

            if h[i % t] == 0: h = op.a_rotate_b(h, B, t)
            else: h = op.a_rotate_b_inv(h, B, t)

        h = op.a_rotate_b_inv(h, k, t)
        self.h = format(en.per_to_int(h, t) % two_exp_s, '0'+str(s * 2)+'x')
        
        return self.h
