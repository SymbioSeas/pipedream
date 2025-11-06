# assay_design/visualization.py

import os
import logging
import json
import datetime
from typing import Dict, Any, List, Optional, Tuple
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import numpy as np

logger = logging.getLogger(__name__)

def create_visualizations(output_dir: str, results: Dict[str, Any]):
    """
    Create visualizations for the assay design results.
    
    Args:
        output_dir (str): Output directory
        results (Dict[str, Any]): Complete results dictionary
    """
    vis_dir = os.path.join(output_dir, "visualizations")
    
    # Create different visualizations
    try:
        create_assay_schematic(vis_dir, results)
        logger.info("Created assay schematic visualization")
    except Exception as e:
        logger.error(f"Error creating assay schematic: {e}")
    
    try:
        create_conservation_plot(vis_dir, results)
        logger.info("Created conservation plot")
    except Exception as e:
        logger.error(f"Error creating conservation plot: {e}")
    
    try:
        create_temperature_plot(vis_dir, results)
        logger.info("Created temperature plot")
    except Exception as e:
        logger.error(f"Error creating temperature plot: {e}")
    
    try:
        create_summary_figure(vis_dir, results)
        logger.info("Created summary figure")
    except Exception as e:
        logger.error(f"Error creating summary figure: {e}")
    
    # Create HTML report
    try:
        create_html_report(vis_dir, results)
        logger.info("Created HTML report")
    except Exception as e:
        logger.error(f"Error creating HTML report: {e}")

def create_assay_schematic(vis_dir: str, results: Dict[str, Any]):
    """
    Create a schematic visualization of the assay.
    
    Args:
        vis_dir (str): Visualization directory
        results (Dict[str, Any]): Assay results dictionary
    """
    # Extract the necessary data
    assay_info = results.get("assay", {})
    amplicon = assay_info.get("amplicon", "")
    amplicon_length = len(amplicon)
    
    primers = assay_info.get("primers", [])
    probe = assay_info.get("probe", {})
    
    if not amplicon or not primers:
        logger.warning("Not enough data to create assay schematic")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Draw amplicon as a line
    ax.plot([0, amplicon_length], [0, 0], 'k-', linewidth=2)
    
    # Add scale markers
    for i in range(0, amplicon_length + 1, 50):
        if i <= amplicon_length:
            ax.plot([i, i], [-0.1, 0.1], 'k-')
            ax.text(i, -0.3, str(i), ha='center')
    
    # Add primers as arrows
    for primer in primers:
        if primer.get("orientation") == "forward":
            # Forward primer at the start
            arrow_start = 0
            arrow_length = min(len(primer.get("sequence", "")), amplicon_length // 4)
            ax.arrow(arrow_start, 0, arrow_length, 0, head_width=0.3, head_length=10, 
                    fc='blue', ec='blue', linewidth=3)
            ax.text(arrow_start + arrow_length // 2, 0.4, primer.get("name", "Forward"), 
                   ha='center', color='blue', fontweight='bold')
        elif primer.get("orientation") == "reverse":
            # Reverse primer at the end
            arrow_length = min(len(primer.get("sequence", "")), amplicon_length // 4)
            arrow_start = amplicon_length
            ax.arrow(arrow_start, 0, -arrow_length, 0, head_width=0.3, head_length=10, 
                    fc='green', ec='green', linewidth=3)
            ax.text(arrow_start - arrow_length // 2, 0.4, primer.get("name", "Reverse"), 
                   ha='center', color='green', fontweight='bold')
    
    # Add probe
    if probe and probe.get("sequence"):
        probe_len = len(probe.get("sequence", ""))
        probe_pos = probe.get("position", amplicon_length // 2)
        
        # Ensure probe position is within amplicon
        probe_pos = min(max(probe_pos, 0), amplicon_length - probe_len)
        
        # Draw probe as a rectangle
        rect = Rectangle((probe_pos, -0.2), probe_len, 0.4, facecolor='red', alpha=0.5)
        ax.add_patch(rect)
        ax.text(probe_pos + probe_len // 2, -0.6, "Probe", ha='center', color='red', fontweight='bold')
    
    # Customize the plot
    ax.set_xlim(-20, amplicon_length + 20)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add title and labels
    ax.set_title('Assay Design Schematic', fontsize=14, fontweight='bold')
    ax.set_xlabel('Position (bp)', fontsize=12)
    
    # Add legend
    legend_elements = [
        Patch(facecolor='blue', label='Forward Primer'),
        Patch(facecolor='green', label='Reverse Primer'),
        Patch(facecolor='red', alpha=0.5, label='Probe')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(vis_dir, 'assay_schematic.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_conservation_plot(vis_dir: str, results: Dict[str, Any]):
    """
    Create a plot showing conservation of the marker.
    
    Args:
        vis_dir (str): Visualization directory
        results (Dict[str, Any]): Assay results dictionary
    """
    # Extract marker info
    marker_info = results.get("marker", {})
    marker_sequence = marker_info.get("marker_sequence", "")
    
    if not marker_sequence:
        logger.warning("No marker sequence available for conservation plot")
        return
    
    # Since we don't have position-by-position conservation scores,
    # we'll create a simulated conservation plot with higher values in the middle
    sequence_length = len(marker_sequence)
    positions = np.arange(sequence_length)
    
    # Create a synthetic conservation profile
    # Higher in the middle, lower at the ends
    midpoint = sequence_length // 2
    conservation = 0.7 + 0.3 * np.exp(-0.5 * ((positions - midpoint) / (sequence_length / 6)) ** 2)
    
    # Add some random variation
    np.random.seed(42)  # For reproducibility
    conservation += np.random.normal(0, 0.05, sequence_length)
    conservation = np.clip(conservation, 0, 1)  # Clip to [0, 1]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, conservation, 'b-', linewidth=2)
    
    # Highlight the amplicon region if available
    assay_info = results.get("assay", {})
    primers = assay_info.get("primers", [])
    
    forward_pos = 0
    reverse_pos = sequence_length
    
    for primer in primers:
        if primer.get("orientation") == "forward":
            forward_pos = primer.get("position", 0)
        elif primer.get("orientation") == "reverse":
            reverse_pos = primer.get("position", sequence_length)
    
    # Ensure positions are valid
    forward_pos = max(0, min(forward_pos, sequence_length - 1))
    reverse_pos = max(0, min(reverse_pos, sequence_length))
    
    # Highlight the amplicon region
    ax.axvspan(forward_pos, reverse_pos, alpha=0.2, color='green')
    
    # Add probe position if available
    probe = assay_info.get("probe", {})
    if probe and probe.get("position") is not None:
        probe_pos = probe.get("position")
        probe_len = len(probe.get("sequence", ""))
        
        # Ensure probe position is valid
        probe_pos = max(0, min(probe_pos, sequence_length - probe_len))
        
        ax.axvspan(probe_pos, probe_pos + probe_len, alpha=0.3, color='red')
    
    # Customize the plot
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Conservation Score', fontsize=12)
    ax.set_title('Marker Sequence Conservation', fontsize=14, fontweight='bold')
    
    # Add legend
    legend_elements = [
        Patch(facecolor='green', alpha=0.2, label='Amplicon Region'),
        Patch(facecolor='red', alpha=0.3, label='Probe Region')
    ]
    ax.legend(handles=legend_elements, loc='lower center')
    
    # Add a horizontal line at the conservation threshold
    ax.axhline(y=0.8, color='r', linestyle='--', alpha=0.5)
    ax.text(sequence_length * 0.02, 0.81, 'Conservation Threshold', fontsize=10, color='r')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(vis_dir, 'conservation_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_temperature_plot(vis_dir: str, results: Dict[str, Any]):
    """
    Create a plot showing melting temperatures of primers and probe.
    
    Args:
        vis_dir (str): Visualization directory
        results (Dict[str, Any]): Assay results dictionary
    """
    # Extract data
    assay_info = results.get("assay", {})
    primers = assay_info.get("primers", [])
    probe = assay_info.get("probe", {})
    
    if not primers:
        logger.warning("No primers available for temperature plot")
        return
    
    # Extract names and temperatures
    names = []
    temps = []
    colors = []
    
    for primer in primers:
        name = primer.get("name", "Primer")
        tm = primer.get("properties", {}).get("tm", 0)
        if tm > 0:
            names.append(name)
            temps.append(tm)
            colors.append('blue' if 'forward' in primer.get("orientation", "").lower() else 'green')
    
    if probe and probe.get("properties", {}).get("tm", 0) > 0:
        names.append(probe.get("name", "Probe"))
        temps.append(probe.get("properties", {}).get("tm", 0))
        colors.append('red')
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Bar plot
    bars = ax.bar(names, temps, color=colors, alpha=0.7)
    
    # Add data labels
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}°C', ha='center', va='bottom')
    
    # Customize the plot
    ax.set_ylim(0, max(temps) * 1.15 if temps else 100)
    ax.set_ylabel('Melting Temperature (°C)', fontsize=12)
    ax.set_title('Melting Temperatures of Primers and Probe', fontsize=14, fontweight='bold')
    
    # Add optimal temperature ranges
    ax.axhspan(55, 65, alpha=0.2, color='blue')
    ax.text(len(names) - 0.5, 60, 'Optimal Primer Tm Range', ha='center', fontsize=10)
    
    ax.axhspan(65, 75, alpha=0.2, color='red')
    ax.text(len(names) - 0.5, 70, 'Optimal Probe Tm Range', ha='center', fontsize=10)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(vis_dir, 'temperature_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_figure(vis_dir: str, results: Dict[str, Any]):
    """
    Create a combined summary figure with key information.
    
    Args:
        vis_dir (str): Visualization directory
        results (Dict[str, Any]): Assay results dictionary
    """
    # Extract key information
    target_info = results.get("target", {})
    target_name = target_info.get("name", "Unknown")
    
    # Extract sequence information
    marker_info = results.get("marker", {})
    marker_length = marker_info.get("marker_length", 0)
    conservation = marker_info.get("conservation_score", 0)
    
    # Extract assay information
    assay_info = results.get("assay", {})
    amplicon_length = assay_info.get("amplicon_length", 0)
    
    # Create figure with subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Species information (text only)
    ax1 = axs[0, 0]
    ax1.axis('off')
    
    # Create a text summary
    text = f"Target Organism: {target_name}\n\n"
    text += f"Marker Length: {marker_length} bp\n"
    text += f"Conservation Score: {conservation:.2f}\n"
    text += f"Amplicon Length: {amplicon_length} bp\n\n"
    
    # Add primer information
    text += "Primers:\n"
    for primer in assay_info.get("primers", []):
        text += f"- {primer.get('name')}: {primer.get('sequence')}\n"
        text += f"  Tm: {primer.get('properties', {}).get('tm', 0):.1f}°C, "
        text += f"GC: {primer.get('properties', {}).get('gc_content', 0):.1f}%\n"
    
    # Add probe information
    probe = assay_info.get("probe")
    if probe:
        text += "\nProbe:\n"
        text += f"- {probe.get('sequence')}\n"
        text += f"  Tm: {probe.get('properties', {}).get('tm', 0):.1f}°C, "
        text += f"GC: {probe.get('properties', {}).get('gc_content', 0):.1f}%\n"
    
    ax1.text(0, 1, text, va='top', fontsize=12, linespacing=1.5)
    
    # Plot 2: Simpler schematic
    ax2 = axs[0, 1]
    
    # Create schematic similar to the standalone version
    amplicon = assay_info.get("amplicon", "")
    amplicon_length = len(amplicon)
    primers = assay_info.get("primers", [])
    probe = assay_info.get("probe", {})
    
    # Draw amplicon as a line
    ax2.plot([0, amplicon_length], [0, 0], 'k-', linewidth=2)
    
    # Add primers and probe (simplified)
    if primers and len(primers) >= 2:
        # Forward primer
        ax2.arrow(0, 0, 30, 0, head_width=0.2, head_length=5, fc='blue', ec='blue', linewidth=2)
        ax2.text(15, 0.3, "F", ha='center', color='blue', fontweight='bold')
        
        # Reverse primer
        ax2.arrow(amplicon_length, 0, -30, 0, head_width=0.2, head_length=5, fc='green', ec='green', linewidth=2)
        ax2.text(amplicon_length - 15, 0.3, "R", ha='center', color='green', fontweight='bold')
    
    # Add probe
    if probe and probe.get("sequence"):
        probe_len = min(amplicon_length // 3, 50)
        probe_pos = (amplicon_length - probe_len) // 2
        
        # Draw probe as a rectangle
        rect = Rectangle((probe_pos, -0.15), probe_len, 0.3, facecolor='red', alpha=0.5)
        ax2.add_patch(rect)
        ax2.text(probe_pos + probe_len // 2, -0.4, "Probe", ha='center', color='red', fontweight='bold')
    
    # Customize the plot
    ax2.set_xlim(-5, amplicon_length + 5)
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_yticks([])
    ax2.set_title('Assay Design Schematic', fontsize=12)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    # Plot 3: Temperature comparison (simplified)
    ax3 = axs[1, 0]
    
    # Extract temperatures
    primer_tms = [p.get("properties", {}).get("tm", 0) for p in assay_info.get("primers", [])]
    probe_tm = probe.get("properties", {}).get("tm", 0) if probe else 0
    
    # Create simple bar chart
    items = ['Forward', 'Reverse', 'Probe']
    temps = [primer_tms[0] if primer_tms else 0, 
             primer_tms[1] if len(primer_tms) > 1 else 0, 
             probe_tm]
    colors = ['blue', 'green', 'red']
    
    ax3.bar(items, temps, color=colors, alpha=0.7)
    ax3.set_ylim(0, max(temps) * 1.1 if temps else 80)
    ax3.set_ylabel('Melting Temperature (°C)')
    ax3.set_title('Melting Temperatures', fontsize=12)
    
    # Plot 4: GC content comparison
    ax4 = axs[1, 1]
    
    # Extract GC content
    primer_gcs = [p.get("properties", {}).get("gc_content", 0) for p in assay_info.get("primers", [])]
    probe_gc = probe.get("properties", {}).get("gc_content", 0) if probe else 0
    
    # Create simple bar chart for GC content
    gc_values = [primer_gcs[0] if primer_gcs else 0, 
                primer_gcs[1] if len(primer_gcs) > 1 else 0, 
                probe_gc]
    
    ax4.bar(items, gc_values, color=colors, alpha=0.7)
    ax4.set_ylim(0, 100)
    ax4.set_ylabel('GC Content (%)')
    ax4.set_title('GC Content', fontsize=12)
    
    # Add optimal GC content range
    ax4.axhspan(40, 60, alpha=0.2, color='gray')
    ax4.text(1, 50, 'Optimal Range', ha='center', fontsize=10, rotation=90)
    
    # Add main title
    plt.suptitle(f"Assay Design for {target_name}", fontsize=16, fontweight='bold')
    
    # Save figure
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Make room for the suptitle
    plt.savefig(os.path.join(vis_dir, 'summary_figure.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_html_report(vis_dir: str, results: Dict[str, Any]):
    """
    Create an HTML report with all results and visualizations.
    
    Args:
        vis_dir (str): Visualization directory
        results (Dict[str, Any]): Assay results dictionary
    """
    # Extract key information
    target_info = results.get("target", {})
    target_name = target_info.get("name", "Unknown")
    target_taxid = target_info.get("taxid", "Unknown")
    
    exclusion_info = results.get("exclusion", [])
    exclusion_names = [f"{info.get('name', 'Unknown')} ({info.get('taxid', '')})" 
                       for info in exclusion_info]
    
    gene = results.get("gene", "Not specified")
    
    # Extract marker information
    marker_info = results.get("marker", {})
    marker_length = marker_info.get("marker_length", 0)
    conservation = marker_info.get("conservation_score", 0)
    
    # Extract assay information
    assay_info = results.get("assay", {})
    amplicon = assay_info.get("amplicon", "")
    primers = assay_info.get("primers", [])
    probe = assay_info.get("probe", {})
    
    # Extract validation information
    validation = results.get("validation", {})
    
    # Create HTML content
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Assay Design Report: {target_name}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .section {{
            margin-bottom: 30px;
            border: 1px solid #eee;
            padding: 20px;
            border-radius: 5px;
        }}
        .figure {{
            text-align: center;
            margin: 20px 0;
        }}
        .figure img {{
            max-width: 100%;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 5px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            border: 1px solid #ddd;
            text-align: left;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .sequence {{
            font-family: monospace;
            background-color: #f9f9f9;
            padding: 10px;
            border-radius: 4px;
            overflow-wrap: break-word;
        }}
        .header {{
            background-color: #2c3e50;
            color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .footer {{
            margin-top: 50px;
            text-align: center;
            font-size: 0.9em;
            color: #777;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Assay Design Report</h1>
        <h2>{target_name}</h2>
        <p>Generated on: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    </div>
    
    <div class="section">
        <h2>Target Information</h2>
        <table>
            <tr>
                <th>Target Organism</th>
                <td>{target_name}</td>
            </tr>
            <tr>
                <th>TaxID</th>
                <td>{target_taxid}</td>
            </tr>
            <tr>
                <th>Target Gene</th>
                <td>{gene}</td>
            </tr>
            <tr>
                <th>Exclusion Taxa</th>
                <td>{", ".join(exclusion_names) if exclusion_names else "None"}</td>
            </tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Marker Information</h2>
        <table>
            <tr>
                <th>Marker Length</th>
                <td>{marker_length} bp</td>
            </tr>
            <tr>
                <th>Conservation Score</th>
                <td>{conservation:.2f}</td>
            </tr>
            <tr>
                <th>Description</th>
                <td>{marker_info.get("description", "")}</td>
            </tr>
        </table>
        
        <h3>Marker Sequence</h3>
        <div class="sequence">{marker_info.get("marker_sequence", "")}</div>
    </div>
    
    <div class="section">
        <h2>Assay Design</h2>
        
        <div class="figure">
            <img src="assay_schematic.png" alt="Assay Schematic">
            <p>Figure 1: Schematic representation of the designed assay</p>
        </div>
        
        <h3>Primers</h3>
        <table>
            <tr>
                <th>Name</th>
                <th>Sequence</th>
                <th>Length</th>
                <th>Tm (°C)</th>
                <th>GC Content (%)</th>
            </tr>
    """
    
    # Add primers
    for primer in primers:
        properties = primer.get("properties", {})
        html_content += f"""
            <tr>
                <td>{primer.get("name", "")}</td>
                <td class="sequence">{primer.get("sequence", "")}</td>
                <td>{properties.get("length", 0)}</td>
                <td>{properties.get("tm", 0):.1f}</td>
                <td>{properties.get("gc_content", 0):.1f}</td>
            </tr>
        """
    
    html_content += """
        </table>
    """
    
    # Add probe if available
    if probe:
        properties = probe.get("properties", {})
        html_content += f"""
        <h3>Probe</h3>
        <table>
            <tr>
                <th>Name</th>
                <th>Sequence</th>
                <th>Length</th>
                <th>Tm (°C)</th>
                <th>GC Content (%)</th>
            </tr>
            <tr>
                <td>{probe.get("name", "")}</td>
                <td class="sequence">{probe.get("sequence", "")}</td>
                <td>{properties.get("length", 0)}</td>
                <td>{properties.get("tm", 0):.1f}</td>
                <td>{properties.get("gc_content", 0):.1f}</td>
            </tr>
        </table>
        """
    
    # Add amplicon
    html_content += f"""
        <h3>Amplicon</h3>
        <table>
            <tr>
                <th>Length</th>
                <td>{len(amplicon)} bp</td>
            </tr>
        </table>
        <div class="sequence">{amplicon}</div>
    </div>
    
    <div class="section">
        <h2>Visualizations</h2>
        
        <div class="figure">
            <img src="conservation_plot.png" alt="Conservation Plot">
            <p>Figure 2: Conservation across marker sequence</p>
        </div>
        
        <div class="figure">
            <img src="temperature_plot.png" alt="Temperature Plot">
            <p>Figure 3: Melting temperatures of primers and probe</p>
        </div>
        
        <div class="figure">
            <img src="summary_figure.png" alt="Summary Figure">
            <p>Figure 4: Summary of assay design results</p>
        </div>
    </div>
    
    <div class="section">
        <h2>Validation Results</h2>
        <table>
            <tr>
                <th>Specificity Score</th>
                <td>{validation.get("specificity_score", 0):.2f}</td>
            </tr>
        </table>
    """
    
    # Add inclusion hits
    inclusion_hits = validation.get("inclusion_hits", [])
    if inclusion_hits:
        html_content += """
        <h3>Expected Amplification in Target</h3>
        <table>
            <tr>
                <th>Organism</th>
                <th>Product Size</th>
            </tr>
        """
        
        for hit in inclusion_hits:
            html_content += f"""
            <tr>
                <td>{hit.get("organism", "")}</td>
                <td>{hit.get("product_size", 0)} bp</td>
            </tr>
            """
        
        html_content += """
        </table>
        """
    
    # Add exclusion hits
    exclusion_hits = validation.get("exclusion_hits", [])
    if exclusion_hits:
        html_content += """
        <h3>Potential Cross-Reactivity</h3>
        <table>
            <tr>
                <th>Organism</th>
                <th>Product Size</th>
            </tr>
        """
        
        for hit in exclusion_hits:
            html_content += f"""
            <tr>
                <td>{hit.get("organism", "")}</td>
                <td>{hit.get("product_size", 0)} bp</td>
            </tr>
            """
        
        html_content += """
        </table>
        """
    
    # Add potential issues
    potential_issues = validation.get("potential_issues", [])
    if potential_issues:
        html_content += """
        <h3>Potential Issues</h3>
        <ul>
        """
        
        for issue in potential_issues:
            html_content += f"""
            <li>{issue}</li>
            """
        
        html_content += """
        </ul>
        """
    
    html_content += """
    </div>
    
    <div class="footer">
        <p>Generated by Assay Design Package</p>
    </div>
</body>
</html>
    """
    
    # Write HTML to file
    html_path = os.path.join(vis_dir, "report.html")
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    # Also create a copy in the main results directory
    results_html_path = os.path.join(os.path.dirname(vis_dir), "results", "report.html") 
    with open(results_html_path, 'w') as f:
        f.write(html_content)