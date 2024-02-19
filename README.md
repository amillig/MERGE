# MERGE

MERGE represents a method that combines direct coupling analysis and machine learning techniques to predict a protein's fitness from sequence. It requires a binary parameter file outputted by [plmc](https://github.com/debbiemarkslab/plmc/tree/master) and sequence-fitness pairs.

![MERGE_git](https://github.com/amillig/MERGE/assets/58852023/f3da6124-5bee-41a2-b4be-a9c8cd0c4947)

<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div><div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<h1 id="Prerequisites">Prerequisites<a class="anchor-link" href="#Prerequisites">&#182;</a></h1><h4 id="Download-latest-version-of-UniRef100">Download latest version of UniRef100<a class="anchor-link" href="#Download-latest-version-of-UniRef100">&#182;</a></h4><ol>
<li>Download the latest version of UniRef100
<code>$ wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz</code></li>
<li>Unzip the file to get a fasta file <code>$ gzip -d uniref100.fasta.gz</code></li>
</ol>
<p>For further information see <a href="https://www.uniprot.org/help/downloads">https://www.uniprot.org/help/downloads</a></p>
<h4 id="Installing-HMMER">Installing HMMER<a class="anchor-link" href="#Installing-HMMER">&#182;</a></h4><ol>
<li>Download tarball <code>$ wget http://eddylab.org/software/hmmer/hmmer.tar.gz</code></li>
<li>Unpack the tarball <code>$ tar zxf hmmer.tar.gz</code></li>
<li>Enter directory <code>$ cd hmmer-3.4</code></li>
<li>Set installaion path <code>$ ./configure --prefix /your/install/path</code></li>
<li>Build HMMER <code>$ make</code></li>
<li>Run self tests <code>$ make check</code></li>
<li>Install programs and man pages <code>$ make install</code></li>
<li>Add executable to PATH <code>$ export PATH="/your/install/path/bin:$PATH"</code></li>
</ol>
<p>For further information see <a href="http://hmmer.org/documentation.html">http://hmmer.org/documentation.html</a></p>
<h4 id="Installing-PLMC">Installing PLMC<a class="anchor-link" href="#Installing-PLMC">&#182;</a></h4><ol>
<li>Clone plmc repository <code>$ git clone https://github.com/debbiemarkslab/plmc.git</code></li>
<li><code>$ cd plmc</code></li>
<li>Build with gcc and OpenMP to enable multicore parallelism <code>$ make all-openmp</code></li>
<li>Add executable to PATH <code>$ export PATH="/your/install/path/plmc/bin:$PATH"</code></li>
</ol>
<p>For further information see <a href="https://github.com/debbiemarkslab/plmc">https://github.com/debbiemarkslab/plmc</a></p>
<h4 id="MERGE">MERGE<a class="anchor-link" href="#MERGE">&#182;</a></h4><ol>
<li><code>$ git clone https://github.com/amillig/MERGE.git</code></li>
<li><code>$ cd MERGE</code></li>
<li><code>$ pip install -r requirements.txt</code> <strong><em>AI Requirements file erstellen</em></strong></li>
</ol>
<h1 id="Usage">Usage<a class="anchor-link" href="#Usage">&#182;</a></h1><p>In the following, the most important steps for model construction are briefly described. For details see example given. To build a model of the fitness landscape of a protein and explore it <em>in silico</em>, the following files are required:</p>
<ul>
<li>protein sequence in fasta format</li>
<li>variant-fitness-pairs in csv format</li>
</ul>
<h3 id="1.-Generate-Multiple-Sequence-Alignment-(MSA)-Using-Jackhmmer">1. Generate Multiple Sequence Alignment (MSA) Using Jackhmmer<a class="anchor-link" href="#1.-Generate-Multiple-Sequence-Alignment-(MSA)-Using-Jackhmmer">&#182;</a></h3><p>To generate a multiple sequence alignment, the target sequence must be provided in fasta format and the inclusion threshold (--incT) must be set. Inclusion threshold: Start with $0.5 \cdot L$, with $L$ being the sequence length. If the MSA contains too few sequences, increasing incT is recommended.</p>
<p>Usage: <code>jackhmmer [-options] &lt;seqfile&gt; &lt;seqdb&gt;</code></p>
<h3 id="2.-Post-process-the-MSA">2. Post-process the MSA<a class="anchor-link" href="#2.-Post-process-the-MSA">&#182;</a></h3><p>In a next step, the MSA is being post-processed by</p>
<ol>
<li>excluding all positions, where the wild type sequence has a gap,</li>
<li>excluding all positions that contain more than 30 % gaps,</li>
<li>excluding all sequences that contain more than 50 % gaps.</li>
</ol>
<p>The script <code>sto2a2m.py</code> can be found in the git repository MERGE in the folder <em>scripts</em>.</p>
<p>Usage: <code>python sto2a2m.py -sto &lt;stoFile&gt;</code></p>
<h3 id="3.-Infer-Parameters-for-Potts-Model-Using-PLMC">3. Infer Parameters for Potts-Model Using PLMC<a class="anchor-link" href="#3.-Infer-Parameters-for-Potts-Model-Using-PLMC">&#182;</a></h3><p>Once the a2m file is generated, we proceed with infering the parameters of the statitstical model.</p>
<p>Usage: <code>plmc [options] alignmentfile</code></p>
<h3 id="4.-Construct-and-Explore-Model-of-Fitness-Landscape-Using-MERGE">4. Construct and Explore Model of Fitness Landscape Using MERGE<a class="anchor-link" href="#4.-Construct-and-Explore-Model-of-Fitness-Landscape-Using-MERGE">&#182;</a></h3><p>Finally, a model of the fitness landscape is generated. See the example for details on how to use <em>merge</em>.</p>

</div>
</div>
</div>

# References
“Combining evolutionary probability and machine learning enables data-driven protein engineering with minimized experimental effort” by Alexander-Maurice Illig, Niklas E. Siedhoff, Mehdi D. Davari*, and Ulrich Schwaneberg*

# Author
MERGE was developed and written by Alexander-Maurice Illig at RWTH Aachen University.

# Credits
MERGE uses binary parameter files that are generated with [plmc](https://github.com/debbiemarkslab/plmc/tree/master) written by [John Ingraham](https://github.com/jingraham).
