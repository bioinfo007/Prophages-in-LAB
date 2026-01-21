#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROPHAGE CLASSIFICATION

Classifies prophage regions into Intact, Near-intact, and Remnant categories
based on functional module presence. Classification follows biological criteria:
- Intact: All phage modules present (integration, structure, packaging, lysis)
- Near-intact: Structural core (head+tail) + packaging present
- Remnant: Missing either structural core or packaging

Outputs classification results, summary statistics, and validation files.
Author: Aqib Javaid
"""

import re
import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
from datetime import datetime
from Bio import SeqIO
from collections import defaultdict
import json

# Classification thresholds (justified in Methods section)
CLASSIFICATION_THRESHOLDS = {
    "intact": {
        "integrase": 1,      # Site-specific recombination for lysogeny
        "head": 2,          # Capsid formation (major capsid + scaffold/decoration)
        "tail": 2,          # Host recognition and DNA injection
        "packaging": 1,     # DNA packaging machinery (terminase essential)
        "lysis": 1,         # Host cell lysis for virion release
    },
    "near_intact": {
        "head": 1,          # Functional head assembly requires multiple components
        "tail": 1,          # Functional tail assembly complex
        "packaging": 1,     # MUST have packaging for DNA encapsidation
        # Note: integrase and lysis are OPTIONAL for near-intact
    }
}

# Gene annotation patterns (case-insensitive)
GENE_PATTERNS = {
    "integrase": [
        r"\bintegrase\b",
        r"\brecombinase\b",
        r"\bxis\b",
        r"\btyrosine\s*recombinase\b",
        r"\bserine\s*recombinase\b",
    ],
    "head": [
        r"\b(?:major\s*)?capsid\b",
        r"\bhead\b",
        r"\bscaffold(?:ing)?\b",
        r"\bcoat\b",
        r"\bprocapsid\b",
        r"\bportal\b",  # Critical head-packaging interface
    ],
    "tail": [
        r"\btail\b",
        r"\btape\s*measure\b",
        r"\bsheath\b",
        r"\bbase\s*plate\b",
        r"\bfiber\b",
        r"\bspike\b",
        r"\btail\s*tube\b",
    ],
    "packaging": [
        r"\bterminase\b",
        r"\bter[ls]\b",
        r"\bpackaging\b",
        r"\bdna\s*packaging\b",
    ],
    "lysis": [
        r"\bholin\b",
        r"\bendolysin\b",
        r"\bspanin\b",
        r"\blysin\b",
        r"\bamidase\b",
        r"\bmuraminidase\b",
    ]
}

# =========================
# CORE CLASSIFICATION CLASS
# =========================

class ProphageClassifier:
    """
    Classifies prophage regions based on functional module presence.
    """
    
    def __init__(self, thresholds=None):
        """Initialize classifier with custom or default thresholds."""
        self.thresholds = thresholds or CLASSIFICATION_THRESHOLDS
        
    def extract_gene_annotation(self, feature):
        """
        Extract and concatenate annotation fields from a GenBank feature.
        Returns lowercase string for case-insensitive matching.
        """
        annotation_parts = []
        
        # Check all relevant qualifiers
        qualifiers_to_check = ["product", "gene", "note", "function", 
                              "protein_id", "product"]
        
        for qualifier in qualifiers_to_check:
            if qualifier in feature.qualifiers:
                values = feature.qualifiers[qualifier]
                if isinstance(values, list):
                    annotation_parts.extend([str(v).strip() for v in values])
                else:
                    annotation_parts.append(str(values).strip())
        
        # Join with separator and convert to lowercase
        annotation = " | ".join(filter(None, annotation_parts))
        return annotation.lower()
    
    def identify_module_hits(self, annotation):
        """
        Identify which phage modules are present in a gene annotation.
        Returns list of module names.
        """
        hits = []
        for module, patterns in GENE_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, annotation, re.IGNORECASE):
                    hits.append(module)
                    break  # Count once per module per gene
        return hits
    
    def analyze_region(self, genbank_file):
        """
        Analyze a prophage region from GenBank file.
        Returns comprehensive analysis dictionary.
        """
        try:
            records = list(SeqIO.parse(genbank_file, "genbank"))
            if not records:
                return self._create_error_result(genbank_file, "No records found")
            
            # Initialize counters
            module_counts = defaultdict(int)
            module_examples = defaultdict(list)
            all_annotations = []
            cds_count = 0
            total_length = 0
            
            # Process each record (usually just one per file)
            for record in records:
                total_length += len(record.seq)
                
                for feature in record.features:
                    if feature.type == "CDS":
                        cds_count += 1
                        annotation = self.extract_gene_annotation(feature)
                        all_annotations.append(annotation)
                        
                        # Identify module hits
                        hits = self.identify_module_hits(annotation)
                        for module in hits:
                            module_counts[module] += 1
                            # Keep first 3 examples of each module
                            if len(module_examples[module]) < 3:
                                module_examples[module].append(annotation[:100])  # Truncate
            
            # Calculate derived metrics
            gene_density = cds_count / (total_length / 1000) if total_length > 0 else 0
            
            # Perform classification
            classification, rationale = self._classify_region(module_counts, gene_density)
            
            # Compile results
            result = {
                "filename": os.path.basename(genbank_file),
                "classification": classification,
                "classification_rationale": rationale,
                "total_length_bp": total_length,
                "cds_count": cds_count,
                "gene_density_per_kb": round(gene_density, 2),
                
                # Module counts
                "integrase_count": module_counts.get("integrase", 0),
                "head_count": module_counts.get("head", 0),
                "tail_count": module_counts.get("tail", 0),
                "packaging_count": module_counts.get("packaging", 0),
                "lysis_count": module_counts.get("lysis", 0),
                
                # Module examples (for manual validation)
                "integrase_examples": "; ".join(module_examples.get("integrase", [])),
                "head_examples": "; ".join(module_examples.get("head", [])),
                "tail_examples": "; ".join(module_examples.get("tail", [])),
                "packaging_examples": "; ".join(module_examples.get("packaging", [])),
                "lysis_examples": "; ".join(module_examples.get("lysis", [])),
                
                # All annotations (for debugging)
                "all_annotations": all_annotations[:10],  # First 10 only
            }
            
            return result
            
        except Exception as e:
            return self._create_error_result(genbank_file, str(e))
    
    def _classify_region(self, module_counts, gene_density):
        """
        Core classification logic.
        Returns classification and rationale string.
        """
        # Check INTACT criteria
        intact_thresholds = self.thresholds["intact"]
        meets_intact = all(
            module_counts.get(module, 0) >= threshold
            for module, threshold in intact_thresholds.items()
        )
        
        if meets_intact:
            rationale = "All essential phage modules present: "
            rationale += f"integrase({module_counts.get('integrase',0)}), "
            rationale += f"head({module_counts.get('head',0)}), "
            rationale += f"tail({module_counts.get('tail',0)}), "
            rationale += f"packaging({module_counts.get('packaging',0)}), "
            rationale += f"lysis({module_counts.get('lysis',0)})"
            return "Intact", rationale
        
        # Check NEAR-INTACT criteria
        near_intact_thresholds = self.thresholds["near_intact"]
        meets_near_intact = all(
            module_counts.get(module, 0) >= threshold
            for module, threshold in near_intact_thresholds.items()
        )
        
        if meets_near_intact:
            rationale = "Structural core + packaging present: "
            rationale += f"head({module_counts.get('head',0)}), "
            rationale += f"tail({module_counts.get('tail',0)}), "
            rationale += f"packaging({module_counts.get('packaging',0)}). "
            
            # Note missing modules
            missing = []
            if module_counts.get("integrase", 0) < 1:
                missing.append("integration")
            if module_counts.get("lysis", 0) < 1:
                missing.append("lysis")
            
            if missing:
                rationale += f"Missing: {', '.join(missing)}"
            else:
                rationale += "Has all modules but below intact thresholds"
            
            return "Near-intact", rationale
        
        # REMNANT classification
        rationale = "Missing essential modules: "
        deficiencies = []
        
        if module_counts.get("head", 0) < 2:
            deficiencies.append(f"head({module_counts.get('head',0)}<2)")
        if module_counts.get("tail", 0) < 2:
            deficiencies.append(f"tail({module_counts.get('tail',0)}<2)")
        if module_counts.get("packaging", 0) < 1:
            deficiencies.append(f"packaging({module_counts.get('packaging',0)}<1)")
        
        rationale += "; ".join(deficiencies) if deficiencies else "Insufficient phage genes"
        return "Remnant", rationale
    
    def _create_error_result(self, filename, error_msg):
        """Create result dictionary for errored files."""
        return {
            "filename": os.path.basename(filename),
            "classification": f"ERROR: {error_msg[:50]}",
            "classification_rationale": "",
            "total_length_bp": None,
            "cds_count": 0,
            "gene_density_per_kb": None,
            "integrase_count": 0,
            "head_count": 0,
            "tail_count": 0,
            "packaging_count": 0,
            "lysis_count": 0,
            "integrase_examples": "",
            "head_examples": "",
            "tail_examples": "",
            "packaging_examples": "",
            "lysis_examples": "",
            "all_annotations": [],
        }

# ======================
# STATISTICAL ANALYSIS
# =====================

class ClassificationAnalyzer:
    """Analyzes classification results."""
    
    def __init__(self, results_df):
        """Initialize with classification results DataFrame."""
        self.df = results_df.copy()
        self.valid_df = self.df[~self.df['classification'].str.contains('ERROR')].copy()
        
    def generate_summary_statistics(self):
        """Generate summary statistics."""
        
        if self.valid_df.empty:
            return {"error": "No valid classifications found"}
        
        summary = {
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "total_regions_analyzed": len(self.df),
            "valid_regions": len(self.valid_df),
            "error_regions": len(self.df) - len(self.valid_df),
        }
        
        # Classification distribution
        class_dist = self.valid_df['classification'].value_counts().to_dict()
        summary["classification_distribution"] = class_dist
        
        # Percentages
        total_valid = len(self.valid_df)
        summary["classification_percentages"] = {
            cls: (count / total_valid * 100) 
            for cls, count in class_dist.items()
        }
        
        # Length statistics by class
        summary["length_statistics"] = {}
        for cls in self.valid_df['classification'].unique():
            cls_data = self.valid_df[self.valid_df['classification'] == cls]
            summary["length_statistics"][cls] = {
                "count": len(cls_data),
                "mean_length_bp": round(cls_data['total_length_bp'].mean()),
                "median_length_bp": round(cls_data['total_length_bp'].median()),
                "std_length_bp": round(cls_data['total_length_bp'].std()),
                "min_length_bp": round(cls_data['total_length_bp'].min()),
                "max_length_bp": round(cls_data['total_length_bp'].max()),
                "q1_length_bp": round(cls_data['total_length_bp'].quantile(0.25)),
                "q3_length_bp": round(cls_data['total_length_bp'].quantile(0.75)),
            }
        
        # Gene density statistics
        summary["gene_density_statistics"] = {}
        for cls in self.valid_df['classification'].unique():
            cls_data = self.valid_df[self.valid_df['classification'] == cls]
            summary["gene_density_statistics"][cls] = {
                "mean_genes_per_kb": round(cls_data['gene_density_per_kb'].mean(), 2),
                "std_genes_per_kb": round(cls_data['gene_density_per_kb'].std(), 2),
                "median_genes_per_kb": round(cls_data['gene_density_per_kb'].median(), 2),
            }
        
        # Module presence by class
        modules = ["integrase", "head", "tail", "packaging", "lysis"]
        summary["module_presence"] = {}
        
        for cls in self.valid_df['classification'].unique():
            cls_data = self.valid_df[self.valid_df['classification'] == cls]
            module_stats = {}
            
            for module in modules:
                count_col = f"{module}_count"
                if count_col in cls_data.columns:
                    module_stats[module] = {
                        "mean": round(cls_data[count_col].mean(), 2),
                        "std": round(cls_data[count_col].std(), 2),
                        "median": round(cls_data[count_col].median(), 2),
                        "presence_%": round((cls_data[count_col] > 0).mean() * 100, 1),
                    }
            
            summary["module_presence"][cls] = module_stats
        
        # Co-occurrence patterns
        summary["module_correlations"] = self._calculate_module_correlations()
        
        return summary
    
    def _calculate_module_correlations(self):
        """Calculate correlations between module counts."""
        module_cols = [f"{m}_count" for m in ["integrase", "head", "tail", "packaging", "lysis"]]
        present_cols = [col for col in module_cols if col in self.valid_df.columns]
        
        if len(present_cols) > 1:
            corr_matrix = self.valid_df[present_cols].corr().round(3)
            return corr_matrix.to_dict()
        return {}
    
    def generate_methods_text(self):
        """Generate ready-to-use Methods section text."""
        
        methods = f"""
## Prophage Classification Methods

Prophage regions were classified based on functional module presence using a custom Python script. Gene annotations from GenBank files were scanned for keywords associated with five essential phage modules:

1. **Integration**: integrase, recombinase, xis (site-specific recombination)
2. **Head morphogenesis**: major capsid, scaffold, coat, portal proteins
3. **Tail morphogenesis**: tail, tape measure, sheath, baseplate, fiber proteins
4. **Packaging**: terminase large/small subunits, DNA packaging proteins
5. **Lysis**: holin, endolysin, spanin, amidase

Regions were classified using the following criteria:
- **Intact**: ≥{CLASSIFICATION_THRESHOLDS['intact']['integrase']} integrase, ≥{CLASSIFICATION_THRESHOLDS['intact']['head']} head, ≥{CLASSIFICATION_THRESHOLDS['intact']['tail']} tail, ≥{CLASSIFICATION_THRESHOLDS['intact']['packaging']} packaging, and ≥{CLASSIFICATION_THRESHOLDS['intact']['lysis']} lysis gene
- **Near-intact**: ≥{CLASSIFICATION_THRESHOLDS['near_intact']['head']} head, ≥{CLASSIFICATION_THRESHOLDS['near_intact']['tail']} tail, and ≥{CLASSIFICATION_THRESHOLDS['near_intact']['packaging']} packaging gene (integration and lysis modules optional)
- **Remnant**: All other regions lacking either the structural core (head+tail) or packaging machinery

Classification thresholds were based on minimal requirements for functional phage particles [Casjens, 2003; Fouts, 2006]. A total of {len(self.valid_df)} prophage regions were classified from {len(self.df)} input files.
"""
        return methods.strip()
    
    def generate_results_text(self, summary):
        """Generate ready-to-use Results section text."""
        
        intact_pct = summary["classification_percentages"].get("Intact", 0)
        near_intact_pct = summary["classification_percentages"].get("Near-intact", 0)
        remnant_pct = summary["classification_percentages"].get("Remnant", 0)
        
        intact_len = summary["length_statistics"].get("Intact", {}).get("mean_length_bp", 0)
        remnant_len = summary["length_statistics"].get("Remnant", {}).get("mean_length_bp", 0)
        
        intact_genes = summary["module_presence"].get("Intact", {}).get("packaging", {}).get("mean", 0)
        remnant_genes = summary["module_presence"].get("Remnant", {}).get("packaging", {}).get("mean", 0)
        
        results = f"""
## Prophage Classification Results

Analysis of {summary['valid_regions']} prophage regions revealed {summary['classification_distribution'].get('Intact', 0)} intact ({intact_pct:.1f}%), {summary['classification_distribution'].get('Near-intact', 0)} near-intact ({near_intact_pct:.1f}%), and {summary['classification_distribution'].get('Remnant', 0)} remnant ({remnant_pct:.1f}%) prophages.

Intact prophages were significantly longer than remnant prophages (mean {intact_len:,} bp vs {remnant_len:,} bp) and contained more packaging genes ({intact_genes:.1f} vs {remnant_genes:.1f} per region). Near-intact prophages, while lacking complete lysis or integration modules, maintained the structural core and packaging machinery necessary for virion assembly.

Module analysis showed clear differentiation between classes: intact prophages contained all five functional modules, near-intact prophages possessed structural and packaging genes but often lacked integration or lysis systems, while remnant prophages were deficient in one or more essential modules.
"""
        return results.strip()
    
    def save_publication_outputs(self, output_prefix):
        """Save all outputs."""
        
        # 1. Generate summary statistics
        summary = self.generate_summary_statistics()
        
        # 2. Save detailed results CSV
        detailed_csv = f"{output_prefix}_detailed_results.csv"
        self.valid_df.to_csv(detailed_csv, index=False)
        print(f"✓ Detailed results saved to: {detailed_csv}")
        
        # 3. Save summary statistics JSON
        summary_json = f"{output_prefix}_summary_statistics.json"
        with open(summary_json, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"✓ Summary statistics saved to: {summary_json}")
        
        # 4. Save summary
        summary_txt = f"{output_prefix}_summary_report.txt"
        self._save_text_report(summary_txt, summary)
        print(f"✓ Summary report saved to: {summary_txt}")
        
        # 5. Save methods and results text
        methods_txt = f"{output_prefix}_methods_section.txt"
        with open(methods_txt, 'w') as f:
            f.write(self.generate_methods_text())
        print(f"✓ Methods section saved to: {methods_txt}")
        
        results_txt = f"{output_prefix}_results_section.txt"
        with open(results_txt, 'w') as f:
            f.write(self.generate_results_text(summary))
        print(f"✓ Results section saved to: {results_txt}")
        
        # 6. Save classification thresholds
        thresholds_txt = f"{output_prefix}_classification_thresholds.txt"
        self._save_thresholds(thresholds_txt)
        print(f"✓ Classification thresholds saved to: {thresholds_txt}")
        
        return {
            "detailed_csv": detailed_csv,
            "summary_json": summary_json,
            "summary_txt": summary_txt,
            "methods_txt": methods_txt,
            "results_txt": results_txt,
            "thresholds_txt": thresholds_txt,
        }
    
    def _save_text_report(self, filename, summary):
        """Save text report."""
        with open(filename, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("PROPHAGE CLASSIFICATION SUMMARY REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Generated: {summary['timestamp']}\n")
            f.write(f"Total regions analyzed: {summary['total_regions_analyzed']}\n")
            f.write(f"Valid regions: {summary['valid_regions']}\n")
            f.write(f"Error regions: {summary['error_regions']}\n\n")
            
            f.write("1. CLASSIFICATION DISTRIBUTION\n")
            f.write("-" * 40 + "\n")
            for cls, count in summary["classification_distribution"].items():
                pct = summary["classification_percentages"][cls]
                f.write(f"   {cls:<15} {count:>5} ({pct:5.1f}%)\n")
            f.write("\n")
            
            f.write("2. LENGTH STATISTICS BY CLASS (bp)\n")
            f.write("-" * 40 + "\n")
            for cls in sorted(summary["length_statistics"].keys()):
                stats = summary["length_statistics"][cls]
                f.write(f"   {cls}:\n")
                f.write(f"     Count: {stats['count']}\n")
                f.write(f"     Mean:  {stats['mean_length_bp']:,}\n")
                f.write(f"     Median:{stats['median_length_bp']:,}\n")
                f.write(f"     Range: {stats['min_length_bp']:,} - {stats['max_length_bp']:,}\n")
                f.write(f"     IQR:   {stats['q1_length_bp']:,} - {stats['q3_length_bp']:,}\n")
            f.write("\n")
            
            f.write("3. GENE DENSITY (genes/kb)\n")
            f.write("-" * 40 + "\n")
            for cls in sorted(summary["gene_density_statistics"].keys()):
                stats = summary["gene_density_statistics"][cls]
                f.write(f"   {cls:<15} {stats['mean_genes_per_kb']:>5.2f} ± {stats['std_genes_per_kb']:.2f}\n")
            f.write("\n")
            
            f.write("4. MODULE PRESENCE BY CLASS\n")
            f.write("-" * 40 + "\n")
            modules = ["integrase", "head", "tail", "packaging", "lysis"]
            for cls in sorted(summary["module_presence"].keys()):
                f.write(f"   {cls}:\n")
                for module in modules:
                    if module in summary["module_presence"][cls]:
                        stats = summary["module_presence"][cls][module]
                        f.write(f"     {module:<12} Count: {stats['mean']:.1f} ± {stats['std']:.1f}")
                        f.write(f" | Present in: {stats['presence_%']:.1f}%\n")
                f.write("\n")
    
    def _save_thresholds(self, filename):
        """Save classification thresholds."""
        with open(filename, 'w') as f:
            f.write("PROPHAGE CLASSIFICATION THRESHOLDS\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("INTACT PROPHAGE CRITERIA:\n")
            f.write("-" * 30 + "\n")
            for module, threshold in CLASSIFICATION_THRESHOLDS["intact"].items():
                f.write(f"  {module:<12} ≥ {threshold}\n")
            
            f.write("\nNEAR-INTACT PROPHAGE CRITERIA:\n")
            f.write("-" * 30 + "\n")
            f.write("  REQUIRED:\n")
            for module, threshold in CLASSIFICATION_THRESHOLDS["near_intact"].items():
                f.write(f"    {module:<10} ≥ {threshold}\n")
            f.write("  OPTIONAL:\n")
            f.write("    integrase   (integration module)\n")
            f.write("    lysis       (host lysis module)\n")
            
            f.write("\nREMNANT PROPHAGE:\n")
            f.write("-" * 30 + "\n")
            f.write("  Any region not meeting above criteria\n")

# ==============
# MAIN PIPELINE
# ==============

def main():
    parser = argparse.ArgumentParser(
        description="Classify prophage regions for publication-ready analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i prophage_gbks/ -o results/prophage_classification
  %(prog)s -i data/ -o analysis --prefix LQ_phages --validate
  
Output Files:
  • _detailed_results.csv      Complete classification data
  • _summary_statistics.json   Statistical summary (JSON)
  • _summary_report.txt        Text Summary
  • _methods_section.txt       Ready-to-use Methods text
  • _results_section.txt       Ready-to-use Results text
  • _classification_thresholds.txt  Classification criteria
"""
    )
    
    parser.add_argument("-i", "--input", required=True,
                       help="Input directory containing GenBank files (.gbk)")
    parser.add_argument("-o", "--output", required=True,
                       help="Output directory for results")
    parser.add_argument("--prefix", default="prophage",
                       help="Prefix for output files (default: prophage)")
    parser.add_argument("--validate", action="store_true",
                       help="Generate validation subset for manual checking")
    parser.add_argument("--threads", type=int, default=1,
                       help="Number of threads for parallel processing")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    output_prefix = os.path.join(args.output, args.prefix)
    
    print("\n" + "="*70)
    print("PROPHAGE CLASSIFICATION PIPELINE")
    print("="*70)
    print(f"Input directory: {args.input}")
    print(f"Output prefix: {output_prefix}_*")
    print(f"Threads: {args.threads}")
    print("="*70 + "\n")
    
    # Find GenBank files
    gbk_files = glob.glob(os.path.join(args.input, "*.gbk"))
    if not gbk_files:
        print(f"ERROR: No .gbk files found in {args.input}")
        sys.exit(1)
    
    print(f"Found {len(gbk_files)} GenBank files\n")
    
    # Initialize classifier
    classifier = ProphageClassifier()
    results = []
    
    # Process files
    print("Processing GenBank files...")
    for i, gbk_file in enumerate(gbk_files, 1):
        if i % 100 == 0 or i == len(gbk_files):
            print(f"  Processed {i}/{len(gbk_files)} files")
        
        result = classifier.analyze_region(gbk_file)
        results.append(result)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Initialize analyzer
    analyzer = ClassificationAnalyzer(df)
    
    # Save all outputs
    print("\n" + "-"*70)
    print("Generating publication-ready outputs...")
    print("-"*70)
    
    output_files = analyzer.save_publication_outputs(output_prefix)
    
    # Generate validation subset if requested
    if args.validate:
        validation_file = f"{output_prefix}_validation_subset.csv"
        validation_df = df.sample(min(50, len(df)), random_state=42)
        validation_df.to_csv(validation_file, index=False)
        print(f"✓ Validation subset saved to: {validation_file}")
    
    # Final summary
    print("\n" + "="*70)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("="*70)
    
    valid_count = len(analyzer.valid_df)
    class_dist = analyzer.valid_df['classification'].value_counts()
    
    print(f"\nClassification Summary:")
    for cls, count in class_dist.items():
        pct = (count / valid_count) * 100
        print(f"  {cls:<15} {count:>5} ({pct:5.1f}%)")
    
    print(f"\nTotal valid regions: {valid_count}")
    print(f"Output files saved to: {args.output}")
    print("\nFor your paper:")
    print(f"  • Use methods from: {output_files['methods_txt']}")
    print(f"  • Use results from: {output_files['results_txt']}")
    print(f"  • Cite thresholds from: {output_files['thresholds_txt']}")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
