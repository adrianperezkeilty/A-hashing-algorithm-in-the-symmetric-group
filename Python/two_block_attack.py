
import finite_field_arithmetic
from hash_sn import *
from hash_aux import *

#################################################################################
# 2-block Collision attack
#################################################################################  

# Attributes: m-> message in bytes (length between 2s+1 and 2^s_2 bits)
#             s-> block size 
#             t-> smallest integer s.t t!>2^s
#             p-> largest prime s.t s < p < t!
#
# Introduce difference in B0 and set ~B1=(h1>>B1)<<~h1 to cancel difference
# h0 determined by the reversed half block of B2 and length of message

class block_attack:

    def __init__(self,m,s,t,p):
        
        self.m = m
        self.s = s 
        self.t = t  
        self.p = p

        # Initialize compression counter, difference vector and first states
        self.dif=[]
        self.compression_counter=0
        self.h0=''
        self.h1=''
        self.h2=''

        # Colliding message
        self.m_tilde=''

    # Pad message
    def pad(self):

        s = self.s
        m = self.m
        t = self.t

        # Half block size and original length of the message in bytes
        s_2 = s//2
        length_orig = len(m)//8

        if length_orig < ((2*s)//8) or length_orig > 2**s_2:
            return 'Error'
        
        s = self.s
        length = len(m)
        pad1 = length%s
        pad2 = m[-min(s_2,length):]
        
        if pad1 == s_2:
            m = m+str(10**(s-1))
        elif pad1 > s_2:
            m = m+str(10**((s-1)-(pad1-s_2)))
        else:
            m = m+str(10**((s_2-1)-pad1))
            
        pad2 = "".join(str(elem) for elem in [pad2[-i] for i in range(1,len(pad2)+1)])
        m = m+pad2.rjust(s_2,'0')

        return m

    # First states h_0,h_1,h_2 assuming message has at least two blocks
    def first_states(self,m_pad,B0_per,B1_per,length_after_pad):

        s = self.s
        m = self.m
        t = self.t

        s_2 = s//2
        a = en.int_to_per((len(m)//8)*2**s_2,t)
        b = en.int_to_per(en.bits_to_int(''.join(m_pad[k] for k in range(length_after_pad-1,length_after_pad-s_2-1,-1)))*2**s_2,t)
        self.h0 = en.a_rotate_b_inv(a,b)

        if self.h0[0]%2: self.h1 = en.a_rotate_b_inv(self.h0,B0_per)
        else: self.h1 = en.a_rotate_b(self.h0,B0_per)

        if self.h1[1]%2: self.h2 = en.a_rotate_b_inv(self.h1,B1_per)
        else: self.h2 = en.a_rotate_b(self.h0,B1_per)

        self.compression_counter += 3

        return self.h0,self.h1,self.h2


    # Block value when seen as a factor for the computation of k
    # If block is zero set it to k
    def block_k(self,k,block):

        if block == 0: return k
        else: return block
        

    # Run a 2-block attack
    def block_attack(self):
        
        s = self.s
        t = self.t
        p = self.p
        m = self.m
        field_obj = finite_field_arithmetic.field(p)

        # Pad message    
        m_pad = block_attack.pad(self)

        if m_pad == 'Error':
            return m_pad

        # Blocks
        B0 = en.bits_to_int(m_pad[:s])
        B1 = en.bits_to_int(m_pad[s:2*s])

        length_after_pad = len(m_pad)
        h0,h1,h2 = block_attack.first_states(self,m_pad, en.int_to_per(B0,t), en.int_to_per(B1,t),length_after_pad)

        # Loop through all blocks B0_tilde for the 2-block attack
        for i in range(2**s):
 
            # First modified block
            B0_tilde = i #randint(0,2**s-1)

            # If ~B0=B0 skip
            if B0_tilde == B0: continue

            B0_k = block_attack.block_k(self,1,B0)
            B1_k = block_attack.block_k(self,B0_k,B1)
            B0_tilde_k = block_attack.block_k(self,1,B0_tilde)
            
            B0_tilde_k_inv = finite_field_arithmetic.field.inverse(field_obj,B0_tilde_k)

            result1 = B0_k*B1_k*B0_tilde_k_inv % p
            
            # Second modified block is valid if B1_tilde < 2^s 
            # when obtained with (>>) and *
            
            if result1 < 2**s:
                
                B0_tilde_per = en.int_to_per(B0_tilde,t)
                if h0[0]%2: h1_tilde = en.a_rotate_b_inv(h0,B0_tilde_per)
                else: h1_tilde = en.a_rotate_b(h0,B0_tilde_per)

                # Compute B1_tilde = (h1>>B1)<<h1_tilde = h2<<h1_tilde
                if h1_tilde[1]%2:
                    B1_tilde_per = en.a_rotate_b(h1_tilde,h2)
                    #h2_tilde = a_rotate_b_inv(h1_tilde,B1_tilde_per)
                else:
                    B1_tilde_per = en.a_rotate_b_inv(h1_tilde,h2)
                    #h2_tilde = a_rotate_b(h1_tilde,B1_tilde_per)

                self.compression_counter += 3

                B1_tilde = en.per_to_int(B1_tilde_per,t)
                result2 = B1_tilde
                
                self.dif.append(abs(result1-result2))
                
                if result1 == result2:

                    # Construct the colliding message
                    self.m_tilde=format(B0_tilde,'0'+str(s)+'b')+format(B1_tilde,'0'+str(s)+'b')+m[2*s:]

                    return 'S'
                
        return 'F'


# Run 2-block attack with a given threshold
# Ex: two_block_attack(16,9,362867,100)
def two_block_attack(s,t,p,threshold):
    
    L,N = [],[]
    count_messages,total_compressions = 0,0

    for i in range(threshold):

        
        # Random message with last half block fixed
        rand_two_first_blocks = format(random.randint(0,2**(2*s)-1),'0'+str(2*s)+'b')
        rand_last_half_block = format(random.randint(0,2**(s//2)-1),'0'+str(s//2)+'b')
        rand_m = rand_two_first_blocks + rand_last_half_block

        print('Attack on ',rand_m[:s],rand_m[s:2*s],rand_m[2*s:])

        # Instance of attack
        attack = block_attack(rand_m,s,t,p)

        # Launch attack
        A=attack.block_attack()

        # Update number of compressions and messages processed
        count_messages += 1
        total_compressions += attack.compression_counter

        if A == 'S':
            
            print(''.ljust(50,'-'))
            print('Success after '+str(count_messages)+' random messages')
            print('M  =',attack.m[:s],attack.m[s:2*s],attack.m[2*s:])            
            print('M~ =',attack.m_tilde[:s],attack.m_tilde[s:2*s],attack.m_tilde[2*s:])
            print('Number of compressions:',total_compressions)

            L.append(count_messages)
            N.append(total_compressions)
            
            total_compressions=count_messages=0

    return L,N

        
