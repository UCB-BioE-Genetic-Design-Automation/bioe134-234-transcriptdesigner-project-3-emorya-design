import os
import traceback
import csv
import time
from statistics import mean
from genedesign.seq_utils.translate import Translate
from genedesign.transcript_designer_v5 import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker

def parse_fasta(fasta_file):
    sequences = {}
    current_gene = None
    current_sequence = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                gene_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split('|')[2].split(' ')[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene:
            sequences[current_gene] = ''.join(current_sequence)
    return sequences

def benchmark_proteome(fasta_file):
    designer = TranscriptDesigner()
    designer.initiate()
    proteome = parse_fasta(fasta_file)
    successful_results = []
    error_results = []
    for gene, protein in proteome.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")
            ignores = set()
            transcript = designer.run(protein, ignores)
            successful_results.append({'gene': gene, 'protein': protein, 'transcript': transcript})
        except Exception as e:
            error_results.append({'gene': gene, 'protein': protein, 'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"})
    return successful_results, error_results

def analyze_errors(error_results):
    error_summary = {}
    with open('error_summary.txt', 'w') as f:
        for error in error_results:
            error_message = error['error'].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")
    return error_summary

def validate_transcripts(successful_results):
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    translator = Translate()
    translator.initiate()
    codon_checker = CodonChecker()
    codon_checker.initiate()
    validation_failures = []
    for result in successful_results:
        cds = ''.join(result['transcript'].codons)
        try:
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")
            original_protein = result['protein']
            translated_protein = translator.run(cds)
            if original_protein != translated_protein:
                raise ValueError(f"Translation mismatch: Original {original_protein}, Translated {translated_protein}")
            if not (cds.startswith(("ATG", "GTG", "TTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")
        except ValueError as e:
            validation_failures.append({'gene': result['gene'], 'protein': result['protein'], 'cds': cds, 'site': f"Translation or completeness error: {str(e)}"})
            continue
        transcript_dna = result['transcript'].rbs.utr.upper() + cds
        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if not passed_hairpin:
            formatted_hairpin = hairpin_string.replace('\n', ' ').replace('"', "'")
            validation_failures.append({'gene': result['gene'], 'protein': result['protein'], 'cds': transcript_dna, 'site': f"Hairpin detected: {formatted_hairpin}"})
        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if not passed_forbidden:
            validation_failures.append({'gene': result['gene'], 'protein': result['protein'], 'cds': transcript_dna, 'site': f"Forbidden sequence: {forbidden_site}"})
        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if not passed_promoter:
            validation_failures.append({'gene': result['gene'], 'protein': result['protein'], 'cds': transcript_dna, 'site': f"Constitutive promoter detected: {found_promoter}" if found_promoter else "Constitutive promoter detected"})
        codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(result['transcript'].codons)
        if not codons_above_board:
            validation_failures.append({'gene': result['gene'], 'protein': result['protein'], 'cds': cds, 'site': f"Codon usage check failed: Diversity={codon_diversity}, Rare Codons={rare_codon_count}, CAI={cai_value}"})
    return validation_failures

def write_validation_report(validation_failures):
    with open('validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures):
    total_validation_failures = len(validation_failures)
    checker_failures = {
        'Forbidden Sequence Checker': 0, 'Hairpin Checker': 0,
        'Codon Usage Checker': 0, 'Promoter Checker': 0,
        'Translation/Completeness Checker': 0
    }
    for failure in validation_failures:
        site = failure['site']
        if "Forbidden sequence" in site: checker_failures['Forbidden Sequence Checker'] += 1
        elif "Hairpin detected" in site: checker_failures['Hairpin Checker'] += 1
        elif "Codon usage check failed" in site: checker_failures['Codon Usage Checker'] += 1
        elif "Constitutive promoter detected" in site: checker_failures['Promoter Checker'] += 1
        elif "Translation or completeness error" in site: checker_failures['Translation/Completeness Checker'] += 1
    with open('summary_report.txt', 'w') as f:
        f.write(f"Total genes processed: {total_genes}\n")
        f.write(f"Parsing runtime: {parsing_time:.2f} seconds\n")
        f.write(f"Execution runtime: {execution_time:.2f} seconds\n")
        f.write(f"Total exceptions: {sum(errors_summary.values())}\n")
        if errors_summary:
            f.write(f"\nTop 3 most common exceptions:\n")
            for error, count in sorted(errors_summary.items(), key=lambda x: x[1], reverse=True)[:3]:
                f.write(f"- {error}: {count} occurrences\n")
        else:
            f.write("No exceptions encountered.\n")
        f.write(f"\nTotal validation failures: {total_validation_failures}\n")
        f.write("\nValidation Failures by Checker:\n")
        for checker, count in checker_failures.items():
            f.write(f"- {checker}: {count} occurrences\n")

def run_benchmark(fasta_file):
    start_time = time.time()
    parsing_start = time.time()
    successful_results, error_results = benchmark_proteome(fasta_file)
    parsing_time = time.time() - parsing_start
    errors_summary = analyze_errors(error_results)
    validation_start = time.time()
    validation_failures = validate_transcripts(successful_results)
    execution_time = time.time() - validation_start
    write_validation_report(validation_failures)
    total_genes = len(successful_results) + len(error_results)
    generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures)

if __name__ == "__main__":
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    run_benchmark(fasta_file)
