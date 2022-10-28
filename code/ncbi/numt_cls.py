def get_seq(row):
	sample_start=row['upstream_size']
	if sample_start<0:
		sample_start=0
	sample_end=sample_start+row['sample_size']
	return row['sequence'][sample_start:sample_end]