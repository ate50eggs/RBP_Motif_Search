import sys
import MySQLdb as mdb
from sys import maxint
from Bio import SeqIO


class MotifMatch:

	def __init__(self, id, score, pos, seq, symb, spec):
		self.motifId = id
		self.motifScore = score
		self.position = pos
		self.sequence = seq
		self.symbol = symb
		self.species = spec

	def __cmp__(self, other):
		if self.motifScore == other.motifScore:
			return 0
		elif self.motifScore > other.motifScore:
			return -1
		else:
			return 1

	def __str__(self):
		return "%s\t%s\t%s\t%s\t%s\t%s" % (self.motifId, self.motifScore,
						self.position, self.sequence, self.symbol, self.species)

        def getData(self):
                return (self.motifId, self.motifScore, self.position, self.sequence, self.symbol, self.species)


class SevenmerPreloader:

	def __init__(self, uname, psswd):
		# query structures

                self.q1 = """\n
                select pwm.motif_id, sum(log(probability))
                from position_weight_matrix pwm 
                where (\n"""
                self.q2 = """
                (position = %i AND base = '%s')"""
                self.q3 = """
                )\ngroup by motif_id\nhaving count(0) = %i\norder by sum(log(probability)) desc;\n"""
                self.con = mdb.connect('127.0.0.1', uname, psswd, 'rbp_motifs', unix_socket='/var/lib/mysql/mysql.sock')

		self.insertCmd = "INSERT INTO precomputed_motif_scores VALUES ('%s', '%s', '%s');"	


	def getNmers(self, i):
		bases = ['A', 'C', 'G', 'U']

		if i == 0:
			return bases
		else:
			prevNmers = self.getNmers(i - 1)
			nmers = []
			for b in bases:
				for n in prevNmers:
					nmers.append(n + b)
			return nmers

	def computeScoresAndLoad(self):

		cursor = self.con.cursor()
		insertion_cursor = self.con.cursor()
		sevenMers = self.getNmers(6)
		motifLen = 7 
		for nmer in sevenMers:
			
       	                posSpecific = []
                       	for n in range(len(nmer)):
                       	        posSpecific.append(self.q2 % (n + 1, nmer[n]))
                       	q2all = " or\n".join(posSpecific)
                       	query = self.q1 + q2all + self.q3 % (motifLen)
			cursor.execute(query)
			res = cursor.fetchone()

			while res is not None:
				(id, score) = res
				insertion_cursor.execute(self.insertCmd % (id, nmer, score))
				self.con.commit()
				res = cursor.fetchone()

class MotifSearch:

	def __init__(self):
	
		# query structures

		self.q = """
		SELECT gi.motif_id, gi.gene_symbol, gi.Species, pms.lod_score
                FROM precomputed_motif_scores pms 
                JOIN gene_info gi on pms.id = gi.motif_id 
                WHERE Species in %s AND pms.sevenmer = '%s';
		"""
		self.validChars = ['A', 'C', 'G', 'T', 'U']
		self.con = mdb.connect('localhost', 'motif_search', '', 'rbp_motifs', unix_socket='/var/lib/mysql/mysql.sock')

	def search(self, nmer, motifLen, hitCount, species):

		# if hitCount == 0, get all of the hits. otherwise get the specified number of best hits

		nmer = nmer.upper()
		for s in nmer:
			if s not in self.validChars:
				return False	

		nmer = nmer.replace('T', 'U')
		cursor = self.con.cursor()		
		hits = []

		posSpecific = []
		query = self.q % (species, nmer)

		cursor.execute(query)
		res = cursor.fetchone()

		while res:
			(id, symbol, species, score) = res
			hits.append(MotifMatch(id, score, 0, nmer, symbol, species))

			if hitCount > 0:
				if score > -maxint and len(hits) == hitCount:
					hits.sort(reverse=True)	
					hits = hits[1:hitCount]
			res = cursor.fetchone()

		hits.sort()
		return hits

	def search_all_sequences(self, distinct_nmers, species, motifLen):
		all_scores = {}
		i = 0
		for nmer in distinct_nmers.keys():
			i += 1
			#if i % 1000 == 0:
			#	print i
			seq_scores = self.search(nmer, motifLen, 0, species)
			all_scores[nmer] = seq_scores

		return all_scores
			
