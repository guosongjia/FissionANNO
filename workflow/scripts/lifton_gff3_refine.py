import sys
import argparse
import logging
from collections import defaultdict, OrderedDict
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
try:
    from Bio import SeqIO
    BIO_AVAILABLE = True
except Exception:
    BIO_AVAILABLE = False


def setup_logging(log_file: Optional[str], log_level: str) -> None:
    level = getattr(logging, log_level.upper(), logging.INFO)
    logger = logging.getLogger()
    logger.setLevel(level)
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fmt = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if log_file:
        fh = logging.FileHandler(log_file + ".log", mode="w")
        fh.setFormatter(fmt)
        logger.addHandler(fh)


@dataclass
class Feature:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: OrderedDict

    def length(self) -> int:
        return self.end - self.start + 1

    def attr_str(self) -> str:
        parts = []
        for k, v in self.attrs.items():
            parts.append(f"{k}={v}")
        return ";".join(parts)

    def to_line(self) -> str:
        cols = [
            self.seqid,
            self.source,
            self.type,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            self.attr_str(),
        ]
        return "\t".join(cols)


def parse_attributes(attr_field: str) -> OrderedDict:
    od = OrderedDict()
    if not attr_field or attr_field == ".":
        return od
    for part in attr_field.split(";"):
        if not part:
            continue
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        od[k] = v
    return od


def read_fasta_sequences(fasta_path: str) -> Dict[str, str]:
    if BIO_AVAILABLE:
        seqs = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            seqs[record.id] = str(record.seq).upper()
        return seqs
    seqs = {}
    with open(fasta_path, "r") as fh:
        curr = None
        buf = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr is not None:
                    seqs[curr] = "".join(buf).upper()
                curr = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if curr is not None:
            seqs[curr] = "".join(buf).upper()
    return seqs


def parse_gff3(gff3_path: str) -> List[Feature]:
    feats: List[Feature] = []
    with open(gff3_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                pass
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
            try:
                start_i = int(start)
                end_i = int(end)
            except Exception:
                pass
                continue
            feat = Feature(
                seqid=seqid,
                source=source,
                type=ftype,
                start=start_i,
                end=end_i,
                score=score,
                strand=strand,
                phase=phase,
                attrs=parse_attributes(attrs),
            )
            feats.append(feat)
    return feats


def filter_gene_attributes(feat: Feature) -> None:
    keep_keys = {"ID", "Parent", "Name", "product"}
    new_attrs = OrderedDict()
    for k, v in feat.attrs.items():
        if k in keep_keys:
            new_attrs[k] = v
    ordered = OrderedDict()
    for k in ["ID", "Parent", "Name", "product"]:
        if k in new_attrs:
            ordered[k] = new_attrs[k]
    feat.attrs = ordered


def filter_mrna_attributes(feat: Feature) -> None:
    keep_keys = {"ID", "Parent", "Name", "product", "valid_orf", "truncated_orf", "last_five_pct"}
    new_attrs = OrderedDict()
    for k, v in feat.attrs.items():
        if k in keep_keys:
            new_attrs[k] = v
    ordered = OrderedDict()
    for k in ["ID", "Parent", "Name", "product", "valid_orf", "truncated_orf", "last_five_pct"]:
        if k in new_attrs:
            ordered[k] = new_attrs[k]
    feat.attrs = ordered


def filter_exon_attributes(feat: Feature) -> None:
    keep_keys = {"ID", "Parent", "Name", "product"}
    new_attrs = OrderedDict()
    for k, v in feat.attrs.items():
        if k in keep_keys:
            new_attrs[k] = v
    ordered = OrderedDict()
    for k in ["ID", "Parent", "Name", "product"]:
        if k in new_attrs:
            ordered[k] = new_attrs[k]
    feat.attrs = ordered


def filter_cds_attributes(feat: Feature) -> None:
    keep_keys = {"ID", "Parent", "Name", "product"}
    new_attrs = OrderedDict()
    for k, v in feat.attrs.items():
        if k in keep_keys:
            new_attrs[k] = v
    ordered = OrderedDict()
    for k in ["ID", "Parent", "Name", "product"]:
        if k in new_attrs:
            ordered[k] = new_attrs[k]
    feat.attrs = ordered


def build_hierarchy(feats: List[Feature]) -> Tuple[Dict[str, Feature], Dict[str, List[Feature]], Dict[str, List[Feature]], Dict[str, List[Feature]]]:
    genes: Dict[str, Feature] = {}
    mrnas_by_gene: Dict[str, List[Feature]] = defaultdict(list)
    exons_by_transcript: Dict[str, List[Feature]] = defaultdict(list)
    cds_by_transcript: Dict[str, List[Feature]] = defaultdict(list)
    for f in feats:
        if f.type == "gene":
            gid = f.attrs.get("ID")
            if gid:
                genes[gid] = f
        elif f.type in ("mRNA", "pseudogenic_transcript"):
            pid = f.attrs.get("Parent")
            tid = f.attrs.get("ID")
            if pid and tid:
                mrnas_by_gene[pid].append(f)
        elif f.type == "exon":
            pid = f.attrs.get("Parent")
            if pid:
                exons_by_transcript[pid].append(f)
        elif f.type == "CDS":
            pid = f.attrs.get("Parent")
            if pid:
                cds_by_transcript[pid].append(f)
    return genes, mrnas_by_gene, exons_by_transcript, cds_by_transcript


def sort_segments(strand: str, segs: List[Feature]) -> List[Feature]:
    if strand == "+":
        return sorted(segs, key=lambda x: (x.start, x.end))
    return sorted(segs, key=lambda x: (x.start, x.end), reverse=True)


def match_exon_cds(transcript_id: str, gene_id: str, exons: List[Feature], cds_list: List[Feature]) -> bool:
    ex_sorted = sorted(exons, key=lambda x: (x.start, x.end, x.score, x.strand))
    cds_sorted = sorted(cds_list, key=lambda x: (x.start, x.end, x.score, x.strand))
    if len(ex_sorted) != len(cds_sorted):
        logging.warning(f"Gene {gene_id} Transcript {transcript_id} exon/CDS count mismatch")
        return False
    ok = True
    for e, c in zip(ex_sorted, cds_sorted):
        if not (e.start == c.start and e.end == c.end and e.score == c.score and e.strand == c.strand):
            c.seqid = e.seqid
            c.source = e.source
            c.start = e.start
            c.end = e.end
            c.score = e.score
            c.strand = e.strand
            ok = False
    return ok


def recompute_cds_phase(cds_list: List[Feature], strand: str) -> None:
    ordered = sort_segments(strand, cds_list)
    total = 0
    for i, c in enumerate(ordered):
        if i == 0:
            c.phase = "0"
        else:
            phase = (3 - (total % 3)) % 3
            c.phase = str(phase)
        total += c.length()


def verify_cds_phase(cds_list: List[Feature], strand: str, transcript_id: str, gene_id: str) -> bool:
    ordered = sort_segments(strand, cds_list)
    total = 0
    ok = True
    for i, c in enumerate(ordered):
        expected = "0" if i == 0 else str((3 - (total % 3)) % 3)
        if c.phase not in (".", expected):
            logging.warning(f"Gene {gene_id} Transcript {transcript_id} CDS phase invalid; expected {expected} got {c.phase}")
            ok = False
        total += c.length()
    return ok


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def extract_coding_sequence(seq_by_contig: Dict[str, str], contig: str, strand: str, exons: List[Feature]) -> str:
    ordered = sorted(exons, key=lambda x: (x.start, x.end))
    parts = []
    for e in ordered:
        s = e.start - 1
        t = e.end
        parts.append(seq_by_contig[contig][s:t])
    cds_seq = "".join(parts)
    if strand == "-":
        cds_seq = revcomp(cds_seq)
    return cds_seq


def has_valid_orf(seq: str) -> bool:
    if len(seq) < 3:
        return False
    if len(seq) % 3 != 0:
        return False
    if not seq.startswith("ATG"):
        return False
    stops = {"TAA", "TAG", "TGA"}
    if seq[-3:] not in stops:
        return False
    for i in range(3, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in stops:
            return False
    return True


def has_near_terminal_stop(seq: str, frac: float = 0.05) -> bool:
    if len(seq) < 6:
        return False
    stops = {"TAA", "TAG", "TGA"}
    first = None
    for i in range(3, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in stops:
            first = i
            break
    if first is None:
        return False
    tail = len(seq) - (first + 3)
    return tail <= int(len(seq) * frac)

def find_first_internal_stop(seq: str) -> Optional[int]:
    stops = {"TAA", "TAG", "TGA"}
    for i in range(3, len(seq) - 3, 3):
        if seq[i:i+3] in stops:
            return i
    return None

def truncate_to_stop(seq_by_contig: Dict[str, str], contig: str, strand: str, exons: List[Feature], frac: float = 0.05, min_bp: int = 0) -> Tuple[List[Feature], bool]:
    cds_seq = extract_coding_sequence(seq_by_contig, contig, strand, exons)
    idx = find_first_internal_stop(cds_seq)
    if idx is None:
        return exons, False
    tail = len(cds_seq) - (idx + 3)
    threshold = int(len(cds_seq) * frac)
    if min_bp and min_bp > threshold:
        threshold = min_bp
    if tail > threshold:
        return exons, False
    ordered = sort_segments(strand, exons)
    remaining = tail
    kept = [seg for seg in ordered]
    j = len(kept) - 1
    while remaining > 0 and j >= 0:
        e = kept[j]
        span = e.length()
        if remaining < span:
            if strand == "+":
                e.end = e.end - remaining
            else:
                e.start = e.start + remaining
            remaining = 0
        else:
            remaining -= span
            kept.pop(j)
        j -= 1
    return kept, True




def build_cds_from_exons(tid: str, exons: List[Feature], prev_cds: Optional[List[Feature]]) -> List[Feature]:
    cds_new: List[Feature] = []
    prev_sorted: List[Feature] = []
    if prev_cds:
        prev_sorted = sort_segments(prev_cds[0].strand if prev_cds else "+", prev_cds)
    for i, e in enumerate(sort_segments(exons[0].strand if exons else "+", exons)):
        cid = None
        if prev_sorted and i < len(prev_sorted):
            cid = prev_sorted[i].attrs.get("ID")
        if not cid:
            eid = e.attrs.get("ID")
            base = eid if eid else f"{tid}_{i+1}"
            cid = f"cds_{base}"
        new_attrs = OrderedDict()
        new_attrs["ID"] = cid
        new_attrs["Parent"] = tid
        cds_new.append(
            Feature(
                seqid=e.seqid,
                source=e.source,
                type="CDS",
                start=e.start,
                end=e.end,
                score=e.score,
                strand=e.strand,
                phase=".",
                attrs=new_attrs,
            )
        )
    return cds_new


def extend_to_downstream_stop(seq_by_contig: Dict[str, str], contig: str, strand: str, exons: List[Feature], transcript_id: str, gene_id: str) -> Tuple[List[Feature], bool, int]:
    ordered = sort_segments(strand, exons)
    cds_seq = extract_coding_sequence(seq_by_contig, contig, strand, ordered)
    added_total = 0
    contig_seq = seq_by_contig.get(contig, "")
    if not contig_seq:
        return ordered, False, added_total
    stops = {"TAA", "TAG", "TGA"}
    if strand == "+":
        last = ordered[-1]
        logging.info(f"Gene {gene_id} Transcript {transcript_id} downstream extension on plus strand restricted to last exon")
        ext_start = last.end
        if ext_start >= len(contig_seq):
            return ordered, False, added_total
        ext_seq = contig_seq[ext_start:len(contig_seq)]
        rem = (len(cds_seq) % 3)
        boundary_needed = (3 - rem) % 3
        # Check boundary codon formed by existing tail + first extension bases
        if boundary_needed > 0 and len(ext_seq) >= boundary_needed:
            boundary_codon = cds_seq[-rem:] + ext_seq[:boundary_needed]
            if boundary_codon in stops:
                last.end = last.end + boundary_needed
                added_total = boundary_needed
                logging.info(f"Gene {gene_id} Transcript {transcript_id} codon integrity prior to stop scan len_mod3={rem} pad_used={boundary_needed}")
                return ordered, True, added_total
        ext_seq2 = ext_seq[boundary_needed:]
        logging.info(f"Gene {gene_id} Transcript {transcript_id} codon integrity prior to stop scan len_mod3={rem} pad_used={boundary_needed}")
        for k in range(0, (len(ext_seq2) // 3)):
            codon = ext_seq2[3*k:3*k+3]
            if codon in stops:
                needed = boundary_needed + 3*(k+1)
                last.end = last.end + needed
                added_total = needed
                return ordered, True, added_total
        return ordered, False, added_total
    else:
        last = ordered[-1]
        logging.info(f"Gene {gene_id} Transcript {transcript_id} downstream extension on minus strand restricted to last exon")
        ext_end = last.start - 1
        if ext_end <= 0:
            return ordered, False, added_total
        ext_seq_genomic = contig_seq[0:ext_end]
        ext_seq = revcomp(ext_seq_genomic)
        rem = (len(cds_seq) % 3)
        boundary_needed = (3 - rem) % 3
        # Check boundary codon formed by existing tail + first extension bases
        if boundary_needed > 0 and len(ext_seq) >= boundary_needed:
            boundary_codon = cds_seq[-rem:] + ext_seq[:boundary_needed]
            if boundary_codon in stops:
                last.start = last.start - boundary_needed
                added_total = boundary_needed
                logging.info(f"Gene {gene_id} Transcript {transcript_id} codon integrity prior to stop scan len_mod3={rem} pad_used={boundary_needed}")
                return ordered, True, added_total
        ext_seq2 = ext_seq[boundary_needed:]
        logging.info(f"Gene {gene_id} Transcript {transcript_id} codon integrity prior to stop scan len_mod3={rem} pad_used={boundary_needed}")
        for k in range(0, (len(ext_seq2) // 3)):
            codon = ext_seq2[3*k:3*k+3]
            if codon in stops:
                needed = boundary_needed + 3*(k+1)
                last.start = last.start - needed
                added_total = needed
                return ordered, True, added_total
        return ordered, False, added_total


def extend_to_upstream_start(seq_by_contig: Dict[str, str], contig: str, strand: str, exons: List[Feature], transcript_id: str, gene_id: str, upstream_window_bp: int = 0) -> Tuple[List[Feature], bool, int]:
    ordered = sort_segments(strand, exons)
    cds_seq = extract_coding_sequence(seq_by_contig, contig, strand, ordered)
    window = len(cds_seq) // 2
    if upstream_window_bp and upstream_window_bp > 0:
        window = min(window, upstream_window_bp)
    if window < 3:
        return ordered, False, 0
    contig_seq = seq_by_contig.get(contig, "")
    if not contig_seq:
        return ordered, False, 0
    first = ordered[0]
    if strand == "+":
        ext_end = first.start - 1
        ext_start = max(0, ext_end - window)
        if ext_start >= ext_end:
            return ordered, False, 0
        ext_seq = contig_seq[ext_start:ext_end]
    else:
        ext_end = min(len(contig_seq), first.end + window)
        ext_start = first.end
        if ext_start >= ext_end:
            return ordered, False, 0
        ext_seq_genomic = contig_seq[ext_start:ext_end]
        ext_seq = revcomp(ext_seq_genomic)
    rem = (len(cds_seq) % 3)
    boundary_needed = (3 - rem) % 3
    # Scan codons aligned to reading frame relative to boundary
    start_idx = len(ext_seq) - (boundary_needed + 3)
    for i in range(start_idx, -1, -3):
        codon = ext_seq[i:i+3]
        if codon == "ATG":
            needed = len(ext_seq) - i
            if strand == "+":
                first.start = first.start - needed
            else:
                first.end = first.end + needed
            logging.info(f"Gene {gene_id} Transcript {transcript_id} upstream extension found ATG and applied; added={needed}")
            return ordered, True, needed
    return ordered, False, 0

def is_contig_end(contig_len: int, gene: Feature) -> bool:
    if gene.start <= 1:
        return True
    if gene.end >= contig_len:
        return True
    return False


def update_mrna_attr(mrna: Feature, key: str, value: str) -> None:
    mrna.attrs[key] = value


def renumber_feature_ids(mrnas_by_gene: Dict[str, List[Feature]], exons_by_transcript: Dict[str, List[Feature]], cds_by_transcript: Dict[str, List[Feature]]) -> None:
    for gid, mrnas in mrnas_by_gene.items():
        for mrna in mrnas:
            tid = mrna.attrs.get("ID")
            exs = exons_by_transcript.get(tid, [])
            cds = cds_by_transcript.get(tid, [])
            if exs:
                ordered_ex = sort_segments(mrna.strand, exs)
                for i, e in enumerate(ordered_ex):
                    e.attrs["ID"] = f"{tid}:exon:{i+1}"
            if cds:
                ordered_cds = sort_segments(mrna.strand, cds)
                for i, c in enumerate(ordered_cds):
                    c.attrs["ID"] = f"{tid}:CDS:{i+1}"

def normalize_pseudogenic_transcripts(genes: Dict[str, Feature], mrnas_by_gene: Dict[str, List[Feature]]) -> None:
    for gid, mrnas in mrnas_by_gene.items():
        for mrna in mrnas:
            if mrna.type == "pseudogenic_transcript":
                mrna.type = "mRNA"

def process(gff3_path: str, fasta_path: str, output_path: str, log_path: Optional[str],
            last_pct_frac: float = 0.05, last_pct_min_bp: int = 0,
            upstream_window_bp: int = 0, final_phase_recompute: bool = True) -> None:
    feats = parse_gff3(gff3_path)
    seqs = read_fasta_sequences(fasta_path)
    genes, mrnas_by_gene, exons_by_transcript, cds_by_transcript = build_hierarchy(feats)
    stats_by_gene: Dict[str, Dict[str, str]] = {}
    for gid, g in genes.items():
        filter_gene_attributes(g)
    for gid, mrnas in mrnas_by_gene.items():
        for mrna in mrnas:
            tid = mrna.attrs.get("ID")
            contig = mrna.seqid
            exons = exons_by_transcript.get(tid, [])
            cds_list = cds_by_transcript.get(tid, [])
            if not exons:
                continue
            match_ok = True
            # Keep original exon structure to restore if gene remains non-coding
            orig_exons = [Feature(seqid=e.seqid, source=e.source, type=e.type, start=e.start, end=e.end, score=e.score, strand=e.strand, phase=e.phase, attrs=OrderedDict(e.attrs)) for e in exons]
            if gid not in stats_by_gene:
                stats_by_gene[gid] = {
                    "cds_exon_match": "N",
                    "up_trigger": "N",
                    "up_len": "NA",
                    "down_trigger": "N",
                    "down_len": "NA",
                    "valid_orf": "N",
                    "pseudogene": "N",
                    "truncated_orf": "N",
                }
            if cds_list:
                match_ok = match_exon_cds(tid, gid, exons, cds_list)
                verify_cds_phase(cds_list, mrna.strand, tid, gid)
                if not match_ok:
                    cds_by_transcript[tid] = build_cds_from_exons(tid, exons, cds_list)
                    recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                stats_by_gene[gid]["cds_exon_match"] = "Y" if match_ok else "N"
            contig_len = len(seqs.get(contig, ""))
            gene_feat = genes.get(gid)
            if contig_len and gene_feat and is_contig_end(contig_len, gene_feat):
                update_mrna_attr(mrna, "valid_orf", "False")
                update_mrna_attr(mrna, "truncated_orf", "True")
                continue
            cds_seq = extract_coding_sequence(seqs, contig, mrna.strand, exons)
            if has_valid_orf(cds_seq):
                update_mrna_attr(mrna, "valid_orf", "True")
                if mrna.type == "pseudogenic_transcript":
                    mrna.type = "mRNA"
                mrna.start = min(e.start for e in exons)
                mrna.end = max(e.end for e in exons)
                continue
            if cds_list:
                cds_seq_orig = extract_coding_sequence(seqs, contig, mrna.strand, cds_list)
                stops = {"TAA", "TAG", "TGA"}
                if len(cds_seq_orig) >= 3 and (not cds_seq_orig.startswith("ATG")) and (cds_seq_orig[-3:] in stops) and (find_first_internal_stop(cds_seq_orig) is None):
                    update_mrna_attr(mrna, "valid_orf", "False")
                    gene_feat = genes.get(gid)
                    if gene_feat:
                        gene_feat.type = "pseudogene"
                    if cds_list:
                        cds_by_transcript[tid] = []
                    continue
            first_stop = find_first_internal_stop(cds_seq)
            if first_stop is not None:
                tail = len(cds_seq) - (first_stop + 3)
                trunc_threshold = int(len(cds_seq) * last_pct_frac)
                if last_pct_min_bp and last_pct_min_bp > trunc_threshold:
                    trunc_threshold = last_pct_min_bp
                if tail <= trunc_threshold:
                    trunc_exons, trunc_changed = truncate_to_stop(seqs, contig, mrna.strand, exons, frac=last_pct_frac, min_bp=last_pct_min_bp)
                    if trunc_changed:
                        exons_by_transcript[tid] = trunc_exons
                        cds_by_transcript[tid] = build_cds_from_exons(tid, trunc_exons, cds_list if cds_list else None)
                        recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                        mrna.start = min(e.start for e in trunc_exons)
                        mrna.end = max(e.end for e in trunc_exons)
                    update_mrna_attr(mrna, "valid_orf", "True")
                    update_mrna_attr(mrna, "last_five_pct", "True")
                    if mrna.type == "pseudogenic_transcript":
                        mrna.type = "mRNA"
                    continue
                else:
                    gene_feat = genes.get(gid)
                    if gene_feat:
                        gene_feat.type = "pseudogene"
                    update_mrna_attr(mrna, "valid_orf", "False")
                    if cds_list:
                        cds_by_transcript[tid] = []
                    continue
            adjusted = False
            contig_seq = seqs.get(contig, "")
            if contig_seq:
                ex_candidate = exons
                cds_seq2 = extract_coding_sequence(seqs, contig, mrna.strand, ex_candidate)
                if not has_valid_orf(cds_seq2):
                    ext_exons, ext_changed, ext_len = extend_to_downstream_stop(seqs, contig, mrna.strand, ex_candidate, tid, gid)
                    if ext_changed:
                        ex_candidate = ext_exons
                        exons_by_transcript[tid] = ex_candidate
                        cds_by_transcript[tid] = build_cds_from_exons(tid, ex_candidate, cds_list if cds_list else None)
                        recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                        mrna.start = min(e.start for e in ex_candidate)
                        mrna.end = max(e.end for e in ex_candidate)
                        stats_by_gene[gid]["down_trigger"] = "Y"
                        stats_by_gene[gid]["down_len"] = str(ext_len)
                    cds_seq3 = extract_coding_sequence(seqs, contig, mrna.strand, ex_candidate)
                    if (not has_valid_orf(cds_seq3)) and (not cds_seq3.startswith("ATG")):
                        up_exons, up_changed, up_len = extend_to_upstream_start(seqs, contig, mrna.strand, ex_candidate, tid, gid, upstream_window_bp=upstream_window_bp)
                        if up_changed:
                            ex_candidate = up_exons
                            exons_by_transcript[tid] = ex_candidate
                            cds_by_transcript[tid] = build_cds_from_exons(tid, ex_candidate, cds_list if cds_list else None)
                            stats_by_gene[gid]["up_trigger"] = "Y"
                            stats_by_gene[gid]["up_len"] = str(up_len)
                        cds_seq4 = extract_coding_sequence(seqs, contig, mrna.strand, ex_candidate)
                        # Final decision after all extension attempts
                        if has_valid_orf(cds_seq4):
                            update_mrna_attr(mrna, "valid_orf", "True")
                            adjusted = True
                            if mrna.type == "pseudogenic_transcript":
                                mrna.type = "mRNA"
                            recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                            mrna.start = min(e.start for e in ex_candidate)
                            mrna.end = max(e.end for e in ex_candidate)
                            gene_feat = genes.get(gid)
                            if gene_feat and gene_feat.type == "pseudogene":
                                gene_feat.type = "gene"
            elif cds_list and not match_ok:
                cds_by_transcript[tid] = build_cds_from_exons(tid, exons, cds_list)
                recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                mrna.start = min(e.start for e in exons)
                mrna.end = max(e.end for e in exons)
            if adjusted:
                continue
            exs_curr = exons_by_transcript.get(tid, []) or exons
            cds_seq_curr = extract_coding_sequence(seqs, contig, mrna.strand, exs_curr)
            if has_valid_orf(cds_seq_curr):
                update_mrna_attr(mrna, "valid_orf", "True")
                if mrna.type == "pseudogenic_transcript":
                    mrna.type = "mRNA"
                gene_feat = genes.get(gid)
                if gene_feat and gene_feat.type == "pseudogene":
                    gene_feat.type = "gene"
                mrna.start = min(e.start for e in exs_curr)
                mrna.end = max(e.end for e in exs_curr)
            else:
                contig_len = len(seqs.get(contig, ""))
                if contig_len and is_contig_end(contig_len, genes.get(gid, mrna)):
                    update_mrna_attr(mrna, "valid_orf", "False")
                    update_mrna_attr(mrna, "truncated_orf", "True")
                else:
                    if has_near_terminal_stop(cds_seq_curr, frac=last_pct_frac):
                        trunc_exons, trunc_changed = truncate_to_stop(seqs, contig, mrna.strand, exs_curr, frac=last_pct_frac, min_bp=last_pct_min_bp)
                        if trunc_changed:
                            exons_by_transcript[tid] = trunc_exons
                            cds_by_transcript[tid] = build_cds_from_exons(tid, trunc_exons, cds_list if cds_list else None)
                            recompute_cds_phase(cds_by_transcript.get(tid, []), mrna.strand)
                            mrna.start = min(e.start for e in trunc_exons)
                            mrna.end = max(e.end for e in trunc_exons)
                        update_mrna_attr(mrna, "valid_orf", "True")
                        update_mrna_attr(mrna, "last_five_pct", "True")
                        if mrna.type == "pseudogenic_transcript":
                            mrna.type = "mRNA"
                        gene_feat = genes.get(gid)
                        if gene_feat and gene_feat.type == "pseudogene":
                            gene_feat.type = "gene"
                    else:
                        first_internal = find_first_internal_stop(cds_seq_curr)
                        if first_internal is not None:
                            gene_feat = genes.get(gid)
                            if gene_feat:
                                gene_feat.type = "pseudogene"
                            update_mrna_attr(mrna, "valid_orf", "False")
                            # Restore original structure for pseudogene output
                            exons_by_transcript[tid] = orig_exons
                            cds_by_transcript[tid] = []
                            mrna.start = min(e.start for e in orig_exons)
                            mrna.end = max(e.end for e in orig_exons)
                        else:
                            stops2 = {"TAA", "TAG", "TGA"}
                            no_internal2 = (find_first_internal_stop(cds_seq_curr) is None)
                            has_terminal2 = (len(cds_seq_curr) >= 3 and cds_seq_curr[-3:] in stops2)
                            if no_internal2 and has_terminal2:
                                update_mrna_attr(mrna, "valid_orf", "False")
                                gene_feat = genes.get(gid)
                                if gene_feat:
                                    gene_feat.type = "pseudogene"
                                exons_by_transcript[tid] = orig_exons
                                cds_by_transcript[tid] = []
                                mrna.start = min(e.start for e in orig_exons)
                                mrna.end = max(e.end for e in orig_exons)
                            else:
                                update_mrna_attr(mrna, "valid_orf", "False")
            # update per-gene stats after processing this transcript
            final_exs = exons_by_transcript.get(tid, []) or exons
            final_seq = extract_coding_sequence(seqs, contig, mrna.strand, final_exs)
            stats_by_gene[gid]["valid_orf"] = "Y" if has_valid_orf(final_seq) else "N"
            gene_feat = genes.get(gid)
            if gene_feat and gene_feat.type == "pseudogene":
                stats_by_gene[gid]["pseudogene"] = "Y"
    for gid, mrnas in mrnas_by_gene.items():
        spans = []
        for mrna in mrnas:
            tid = mrna.attrs.get("ID")
            exs = exons_by_transcript.get(tid, [])
            if exs:
                spans.append((min(e.start for e in exs), max(e.end for e in exs)))
            else:
                spans.append((mrna.start, mrna.end))
        if spans and gid in genes:
            contig_len = len(seqs.get(genes[gid].seqid, ""))
            if contig_len and not is_contig_end(contig_len, genes[gid]):
                genes[gid].start = min(s for s, _ in spans)
                genes[gid].end = max(e for _, e in spans)

    normalize_pseudogenic_transcripts(genes, mrnas_by_gene)
    renumber_feature_ids(mrnas_by_gene, exons_by_transcript, cds_by_transcript)

    if final_phase_recompute:
        for gid, mrnas in mrnas_by_gene.items():
            for mrna in mrnas:
                tid = mrna.attrs.get("ID")
                cds_list = cds_by_transcript.get(tid, [])
                if cds_list:
                    recompute_cds_phase(cds_list, mrna.strand)

    out_lines: List[str] = ["##gff-version 3"]
    for gid, gene in genes.items():
        out_lines.append(gene.to_line())
        mrnas = mrnas_by_gene.get(gid, [])
        for mrna in mrnas:
            filter_mrna_attributes(mrna)
            out_lines.append(mrna.to_line())
            tid = mrna.attrs.get("ID")
            exs = sort_segments(mrna.strand, exons_by_transcript.get(tid, []))
            for e in exs:
                filter_exon_attributes(e)
                out_lines.append(e.to_line())
            cds_segs = sort_segments(mrna.strand, cds_by_transcript.get(tid, []))
            for c in cds_segs:
                filter_cds_attributes(c)
                out_lines.append(c.to_line())
    with open(output_path, "w") as out:
        out.write("\n".join(out_lines) + "\n")
    if log_path:
        header = "\t".join([
            "gene_id", "cds_exon_match", "up_trigger", "up_len",
            "down_trigger", "down_len", "valid_orf", "pseudogene", "truncated_orf",
        ])
        lines = [header]
        for gid in genes.keys():
            s = stats_by_gene.get(gid, {
                "cds_exon_match": "N",
                "up_trigger": "N",
                "up_len": "NA",
                "down_trigger": "N",
                "down_len": "NA",
                "valid_orf": "N",
                "pseudogene": "N",
                "truncated_orf": "N",
            })
            mrnas = mrnas_by_gene.get(gid, [])
            has_valid_attr = any(m.attrs.get("valid_orf") == "True" for m in mrnas)
            s["valid_orf"] = "Y" if has_valid_attr else "N"
            gene_feat = genes.get(gid)
            s["pseudogene"] = "Y" if (gene_feat and gene_feat.type == "pseudogene") else "N"
            has_trunc_attr = any(m.attrs.get("truncated_orf") == "True" for m in mrnas)
            s["truncated_orf"] = "Y" if has_trunc_attr else "N"
            cols = [
                gid,
                s["cds_exon_match"],
                s["up_trigger"],
                s["up_len"],
                s["down_trigger"],
                s["down_len"],
                s["valid_orf"],
                s["pseudogene"],
                s["truncated_orf"],
            ]
            lines.append("\t".join(cols))
        with open(log_path, "w") as lf:
            lf.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff3", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file", default=None, help="Path to write per-gene stat table")
    parser.add_argument("--log-level", default="INFO", help="Logging level: INFO, WARNING, ERROR")
    parser.add_argument("--last-pct-frac", type=float, default=0.05,
                        help="Tail fraction (of CDS length) within which an internal stop is considered near-terminal")
    parser.add_argument("--last-pct-min-bp", type=int, default=30,
                        help="Absolute floor (bp) added to --last-pct-frac threshold; default 30 protects short genes from over-strict truncation")
    parser.add_argument("--upstream-window-bp", type=int, default=300,
                        help="Cap on upstream search window for ATG; default 300 prevents extension into adjacent upstream genes (pass 0 to disable cap)")
    parser.add_argument("--no-final-phase-recompute", action="store_true",
                        help="Disable C3: final pass that recomputes phase on all CDS")
    args = parser.parse_args()
    setup_logging(args.log_file, args.log_level)
    process(
        args.gff3, args.fasta, args.output, args.log_file,
        last_pct_frac=args.last_pct_frac,
        last_pct_min_bp=args.last_pct_min_bp,
        upstream_window_bp=args.upstream_window_bp,
        final_phase_recompute=not args.no_final_phase_recompute,
    )


if __name__ == "__main__":
    main()