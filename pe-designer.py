import argparse, time, sys, re
from needle import needle

def rev_comp(s): return s.translate(s.maketrans('AATGCRYSWKMBDHV', 'TACGYRWSMKVHDB'))[::-1]

class PeDesigner:
	
	def __init__(self, args):
		ref_file = args.ref_seq
		ind_file = args.ind_seq
		sref = open(ref_file).readlines()
		sind = open(ind_file).readlines()
		if sref[0] != '>' or sind[0] !='>':
			print("[Error] the file is not FASTA format !!")
			sys.exit()
		self.ref_seq = ''.join(sref[1:]).replace('\n','').upper()
		self.ind_seq = ''.join(sind[1:]).replace('\n','').upper()
		self.pam = args.pam.upper()
		self.ref_dir = args.ref_dir
		self.rna_len = args.l
		self.mm_num = args.m
		self.pbs_min = args.pbs_min
		self.pbs_max = args.pbs_max
		self.rtt_min = args.rtt_min
		self.rtt_max = args.rtt_max
		self.nick_min = args.nick_min
		self.nick_max = args.nick_max
		self.seq_list = [] #d -> [seq, info_l, pbs_d, rtt_d, nick_d]
		self.off_d = {}

		self.offinder_input = [self.ref_dir, 'N'*self.rna_len + pam] 

	def set_mutation(self):
		needle_res = needle(self.ref_seq, self.ind_seq, 10, 10, 0.5, 0.5)
		self.mutation_st = -1
		for n, sym in enumerate(needle_res[1]):
			if sym != '|':
				self.mutation_st = n - 1
		self.mutation_ed = -1
		for n, sym in enumerate(needle_res[1][::-1]):
			if sym != '|':
				self.mutation_end = len(needle_res[1]) - n
		if self.mutation_st == -1 or self.mutation_ed == -1:
			print('reference sequence and induced sequence is same!!!')
			sys.exit()
		self.mutation_len = self.mutation_ed - self.mutation_st - 1
		self.needle_res = needle_res

	def find_target(self):
		pattern = 'N'*self.rna_len + self.pam
		pattern = '(' + pattern + ')|(' + rev_comp(pattern) + ')'
		pattern = pattern.replace('N', '[AGTC]').replace('R', '[AG]').replace('W', '[AT]').replace('M', '[AC]').replace('Y', '[CT]').replace('S', '[GC]').replace('K', '[GT]').replace('B', '[CGT]').replace('D', '[AGT]').replace('H', '[ACT]').replace('V', '[ACG]')
		pattern = re.compile(pattern)
		pos = 0
		m = pattern.search(self.ref_seq, pos)
		targets = [['1', []]]
		while m != None:
			seq = m.group()
			self.seq_list.append(SeqInfo(seq))
			p = m.start()
			if mgroup(1) != None:
				strand = '+'
				rel_pos = round((p+(self.rna_len-3))*100/len(self.ref_seq),2)
			else:
				seq = rev_comp(seq)
				strand = '-'
				rel_pos = round((p+(3+len(self.pam)))*100/len(self.ref_seq),2)
			if seq[:self.rna_len] not in self.offinder_input:
				self.offinder_input.append(seq[:self.rna_len])
			gc = round(seq[:self.rna_len].count('C') + seq[:self.rna_len].count('
			targets[-1][1].append([seq, p + 1, rel_pos, strand, ]])
			pos += 1
			m = pattern.search(ref_seq, pos)
			
			

def parse_args():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("ref_seq", type=str, help="Path of reference sequence file in FASTA format")
	parser.add_argument("ind_seq", type=str, help="Path of reference sequence file in FASTA format")
	parser.add_argument("pam", type=str, help="PAM sequence")
	parser.add_argument("ref_dir", type=str, help="Directory of reference genome file, in FASTA or 2bit format")
	parser.add_argument("-l", type=int, help="length of target without PAM", default=20)
	parser.add_argument("-m", type=int, help="Mismatch number", default=2)
	parser.add_argument("--pbs_min", type=int, help="Minimum of PBS length", default = 12)
	parser.add_argument("--pbs_max", type=int, help="Maximum of PBS length", default = 14)
	parser.add_argument("--rtt_min", type=int, help="Minimum of RTT length", default = 10)
	parser.add_argument("--rtt_max", type=int, help="Maximum of RTT length", default = 20)
	parser.add_argument("--nick_min", type=int, help="Minimum of nicking distance", default = 0)
	parser.add_argument("--nick_max", type=int, help="Maximum of nicking distance", default = 100)
	
	return parser.parse_args()


def find_target():
	

def main():
	
	tst = time.time()
	
	ref_file = args.ref_seq
		
	target_list = []

if __name__ == "__main__":
	main()
