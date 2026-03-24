import sys
import os
import csv
import argparse
import random
import numpy as np
from collections import Counter
import hashlib
import math

class BloomFilter:
    def __init__(self, capacity: int, error_rate: float = 0.001):
        n = max(1, int(capacity))
        p = max(1e-9, min(0.2, float(error_rate)))
        ln2 = math.log(2.0)
        m = int(-(n * math.log(p)) / (ln2 * ln2))
        k = max(1, int((m / n) * ln2))
        self.m = max(8, m)
        self.k = k
        self.bytes = bytearray((self.m + 7) // 8)
        self.n = 0

    def _hashes(self, s: str):
        b = s.encode('utf-8', errors='ignore')
        h1 = int.from_bytes(hashlib.md5(b).digest(), 'big')
        h2 = int.from_bytes(hashlib.sha1(b).digest(), 'big')
        for i in range(self.k):
            yield (h1 + i * h2) % self.m

    def add(self, s: str):
        for idx in self._hashes(s):
            self.bytes[idx >> 3] |= (1 << (idx & 7))
        self.n += 1

    def __contains__(self, s: str) -> bool:
        for idx in self._hashes(s):
            if not (self.bytes[idx >> 3] & (1 << (idx & 7))):
                return False
        return True

def get_max_tries(seq_array):
    freq = {}
    for letter in seq_array:
        freq[letter] = freq.get(letter, 0) + 1
    
    n_fact = math.factorial(len(seq_array))
    product = 1
    for count in freq.values():
        product *= math.factorial(count)
    
    return n_fact / product

def generate_unique_random_peptide(peptide, avoid_contains, avoid_add):
    seq_array = np.array(list(peptide))
    max_tries = get_max_tries(seq_array)
    
    tried_in_this_run = set()
    
    for _ in range(int(max_tries) + 1):
        np.random.shuffle(seq_array)
        shuffled_seq = ''.join(seq_array)
        
        # Check if the sequence is globally unique and not tried in this run
        if shuffled_seq not in tried_in_this_run:
            if not avoid_contains(shuffled_seq):
                avoid_add(shuffled_seq)
                return shuffled_seq
            tried_in_this_run.add(shuffled_seq) # Add to tried_in_this_run even if it was a Bloom filter hit
    
    return None

def iter_csv_rows(path: str):
    with open(path, newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            yield row

def count_rows(path: str) -> int:
    with open(path, 'r', encoding='utf-8') as f:
        return max(0, sum(1 for _ in f) - 1)

def build_target_bloom(path: str, capacity: int, logfh, error_rate: float = 0.001):
    seq_bf = BloomFilter(capacity=capacity, error_rate=error_rate)
    pep_bf = BloomFilter(capacity=capacity, error_rate=error_rate)
    added = 0
    #This function will create two bloom filters: the peptide bf, with only unique peptides and skip duplicates
    #and the seq bf, with only those seqs whose peptide was previuosly added to the peptide bf, therefore, being duplicate seqs with
    #different names in the target DB
    for row in iter_csv_rows(path):
        seq = (row.get('seq_name') or '').strip()
        pep = (row.get('peptide') or '').strip()
        if pep in pep_bf:
            seq_bf.add(seq)
            continue
        pep_bf.add(pep)
        added += 1
        if added % 5000000 == 0:
            (logfh.write if logfh else sys.stderr.write)(f"[INFO] Bloom add: {added:,}\n")
    (logfh.write if logfh else sys.stderr.write)(f"[INFO] Bloom construido con {added:,} péptidos\n")
    return seq_bf, pep_bf

def _process_batch(batch, avoid_contains, avoid_add, writer, logfh):
    total = ok = skipped = 0
    for seq, pep in batch:
        total += 1
        dec = generate_unique_random_peptide(
            pep, avoid_contains=avoid_contains, avoid_add=avoid_add,
        )
        if dec is None:
            skipped += 1
            writer.writerow({'seq_name': seq, 'peptide': pep, 'random_peptide': ''})
            (logfh.write if logfh else sys.stderr.write)(f"[WARN] No decoy para {seq} ('{pep}')\n")
            continue
        else:
            writer.writerow({'seq_name': seq, 'peptide': pep, 'random_peptide': dec})
            ok += 1
    (logfh.write if logfh else sys.stderr.write)(f"[INFO] Progreso: total={total:,} ok={ok:,} skipped={skipped:,}\n")
    return total, ok, skipped

def randomize_file(
    input_csv: str,
    output_csv: str,
    seed: int,
    chunksize: int,
    logfh,
):
    np.random.seed(seed)
    random.seed(seed)

    out_fh = open(output_csv, 'w', newline='', encoding='utf-8')
    writer = csv.DictWriter(out_fh, fieldnames=['seq_name', 'peptide', 'random_peptide'])
    writer.writeheader()

    total = ok = skipped = dups = 0
    #total_rows i.e. the bf capacity needs to be multiplied by 2 to accomodate for the maximun number of decoys to add to the bf
    total_rows = (count_rows(input_csv))*2
    (logfh.write if logfh else sys.stderr.write)(f"[INFO] Filas estimadas (sin header): {total_rows:,}\n")
    seq_bf, pep_bf = build_target_bloom(input_csv, capacity=total_rows, logfh=logfh, error_rate=0.001)

    def avoid_contains(x: str) -> bool:
        return (x in pep_bf)

    def avoid_add(x: str):
        pep_bf.add(x)

    batch = []
    for row in iter_csv_rows(input_csv):
        seq = (row.get('seq_name') or '').strip()
        pep = (row.get('peptide') or '').strip()
        if seq in seq_bf:
            total += 1
            dups += 1
            (logfh.write if logfh else sys.stderr.write)(f"[WARN] {seq} ('{pep}') was duplicated in the original target DB\n")
            continue
        batch.append((seq, pep))
        if len(batch) >= chunksize:
            t, o, s = _process_batch(batch, avoid_contains, avoid_add, writer, logfh)
            total += t; ok += o; skipped += s
            batch = []
    if batch:
        t, o, s = _process_batch(batch, avoid_contains, avoid_add, writer, logfh)
        total += t; ok += o; skipped += s

    out_fh.close()
    (logfh.write if logfh else sys.stderr.write)(
        f"[DONE] total={total:,} | ok={ok:,} | skipped={skipped:,} | duplicated_targets={dups:,}\n[OUT] {output_csv}\n"
    )

def main():
    ap = argparse.ArgumentParser(description="Randomiza péptidos in-silico preservando longitud y composición (streaming/chunks)\nEvita colisiones con targets via Bloom filter.\n")
    ap.add_argument("--input", "-i", required=True, help="CSV de entrada con columnas 'seq_name,peptide'")
    ap.add_argument("--output", "-o", default="randomized_python_output.csv", help="CSV de salida")
    ap.add_argument("--log", default=None, help="Archivo de log (opcional)")
    ap.add_argument("--seed", type=int, default=12345, help="Semilla RNG para reproducibilidad")
    ap.add_argument("--chunksize", type=int, default=1000000, help="Tamaño de lote para procesar filas (memoria/IO)")
    args = ap.parse_args()
    
    logfh = open(args.log, 'w', encoding='utf-8') if args.log else None
    
    try:
        randomize_file(
            input_csv=args.input,
            output_csv=args.output,
            seed=args.seed,
            chunksize=max(1, int(args.chunksize)),
            logfh=logfh,
        )
    finally:
        if logfh:
            logfh.close()

if __name__ == '__main__':
    main()
