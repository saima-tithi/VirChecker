import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.lang.ProcessBuilder.Redirect;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class VirChecker {
    private static String read1 = "";
    private static String read2 = "";   
    private static String outputDir = "";
    private static String scaffoldFile = "";
    private static double avgReadLen = 0;
    private static String spadesKmerlen = "default";
    private static int minOverlapCircular = 5000;
    private static int minIdentityCircular = 95;
    private static double salmonReadFraction = 0;
    private static int minSuspiciousLen = 1000;
    
    private static void parseArguments(String[] args) {
        if (args.length == 0) {
            System.exit(1);
        } else {
            for (int i = 0; i < args.length; i++) {
                if (args[i].startsWith("-")) {
                    if ((i + 1) >= args.length) {
                        System.out.println("Missing argument after " + args[i] + ".");
                        System.exit(1);
                    } else {
                        if (args[i].equals("-o")) {
                            outputDir = args[i + 1];
                        } else if (args[i].equals("-1")) {
                            read1 = args[i + 1];
                        } else if (args[i].equals("-2")) {
                            read2 = args[i + 1];
                        } else if (args[i].equals("-scaffold")) {
                            scaffoldFile = args[i + 1];
                        } else if (args[i].equals("-spadeskmer")) {
                            spadesKmerlen = args[i + 1];
                        } else if (args[i].equals("-minOverlapCircular")) {
                            minOverlapCircular = Integer.parseInt(args[i + 1]);
                        } else if (args[i].equals("-minIdentityCircular")) {
                            minIdentityCircular = Integer.parseInt(args[i + 1]);
                        } else if (args[i].equals("-readFrac")) {
                            salmonReadFraction = Double.parseDouble(args[i + 1]);
                        } else if (args[i].equals("-minSuspiciousLen")) {
                            minSuspiciousLen = Integer.parseInt(args[i + 1]);
                        } else {
                            System.out.println("Invalid argument.");
                            System.exit(1);
                        }
                    }
                }
            }
        } // finish parsing arguments
    }
    
    private static void getReadLen() {
        String cmd = "";
        try {          
            if (read1.endsWith(".gz")) {
                cmd = "gzip -dc " + read1 
                        + " | awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}'";
            } else {
                cmd = "awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}' "
                        + read1;
            }
            
            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String str = "";
            while ((str = reader.readLine()) != null) {
                avgReadLen = Double.parseDouble(str);
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        if (avgReadLen == 0) {
            System.out.println("Could not extract average read length from read file.");
            System.exit(1);
        }
        System.out.println("Estimated average read length: " + avgReadLen);
    }
    
    private static void createBed() {
        //move contig.fasta and contig.fasta.fai from spades-res folder to outputDir
        File scaffoldFile = new File(outputDir + "/scaffold-truncated/tmp/spades-res/scaffold.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            String cmd = "";
            try {
                cmd = "cd " + outputDir + "/scaffold-truncated\n"
                        + "rm scaffold.fasta*\n"
                        + "cd tmp/spades-res\n" 
                        + "mv scaffold.fasta " + outputDir + "/scaffold-truncated\n"
                        + "mv scaffold.fasta.fai " + outputDir + "/scaffold-truncated\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/scaffold-truncated/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/scaffold-truncated/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }           
        }       
        
        BufferedReader br = null;
        BufferedWriter bw = null;
        BufferedWriter bwOutLog = null;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/scaffold-truncated/scaffold.fasta.fai"));
            bw = new BufferedWriter(new FileWriter(outputDir + "/scaffold-truncated/scaffold-start-end.bed"));
            bwOutLog = new BufferedWriter(new FileWriter(outputDir + "/output-log.txt", true));
            String str = br.readLine();
            String[] results = str.split("\t");
            int scaffoldLength = Integer.parseInt(results[1].trim());
            
            if (scaffoldLength >= avgReadLen * 2) {
                String scaffoldId = results[0].trim();
                bw.write(scaffoldId + "\t" + 0 + "\t" + (int) Math.ceil(avgReadLen * 1.5) + "\n");
                bw.write(scaffoldId + "\t" + (int) Math.ceil(scaffoldLength - avgReadLen * 1.5) + "\t" + scaffoldLength
                        + "\n");
                bwOutLog.write("Trying to grow scaffold " + scaffoldId + " with length "
                        + scaffoldLength + "\n");
                if (scaffoldLength > 300000) {
                    bwOutLog.write("Length of " + scaffoldId + " is already greater than 300kbp, so stop extending this one.\n");
                }
            }
            br.close();
            bw.close();
            bwOutLog.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void runAlignment() {
        String cmd = "";
        try {
            if (read2.isEmpty()) {
                cmd = "cd " + outputDir + "/scaffold-truncated\n"
                    + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n"
                    + "rm -r salmon-index\n"
                    + "rm -r salmon-res\n"
                    + "rm salmon-mapped.sam\n"
                    + "salmon index -t scaffold-start-end.fasta -i salmon-index\n"
                    + "/usr/bin/time -f \"\t%E Elasped Real Time\" salmon quant -i salmon-index -l A "
                    + "-r " + read1 + " -o salmon-res --writeMappings -p 16 --quasiCoverage "
                    + salmonReadFraction
                    + " | samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n";
            }
            else {
                    cmd = "cd " + outputDir + "/scaffold-truncated\n"
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n"
                        + "rm -r salmon-index\n"
                        + "rm -r salmon-res\n"
                        + "rm salmon-mapped.sam\n"
                        + "salmon index -t scaffold-start-end.fasta -i salmon-index\n"
                        + "/usr/bin/time -f \"\t%E Elasped Real Time\" salmon quant -i salmon-index -l A "
                        + "-1 " + read1 + " -2 " + read2 + " -o salmon-res --writeMappings -p 16 --quasiCoverage "
                        + salmonReadFraction
                        + "| samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n";
            }
            
            FileWriter shellFileWriter = new FileWriter(outputDir + "/scaffold-truncated/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/scaffold-truncated/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log-alignment.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private static void getMappedReads() {     
        String cmd = "";
        try {
            if (read2.isEmpty()) {
                cmd = "cd " + outputDir + "/scaffold-truncated\n" 
                    + "rm -r tmp\n" 
                    + "mkdir tmp\n"
                    + "bash filterbyname.sh in=" + read1
                    + " out=tmp/mapped_reads_1.fastq names="
                    + "salmon-mapped.sam include=t\n";
            }
            else {
                cmd = "cd " + outputDir + "/scaffold-truncated\n" 
                    + "rm -r tmp\n" 
                    + "mkdir tmp\n"
                    + "bash filterbyname.sh in=" + read1 + " in2=" + read2
                    + " out=tmp/mapped_reads_1.fastq out2=tmp/mapped_reads_2.fastq names="
                    + "salmon-mapped.sam include=t\n";
            }           

            FileWriter shellFileWriter = new FileWriter(outputDir + "/scaffold-truncated/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/scaffold-truncated/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private static void runSpades() {
        String cmd = "";
        try {
            if (read2.isEmpty()) {
                if (spadesKmerlen.equals("default")) {
                    cmd = "cd " + outputDir + "/scaffold-truncated/tmp\n" 
                        + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                        + "spades-res -s mapped_reads_1.fastq "
                        + "--trusted-contigs ../scaffold.fasta"
                        + " --only-assembler\n"
                        + "cd spades-res\n"
                        + "samtools faidx scaffolds.fasta\n";
                }
                else {
                    cmd = "cd " + outputDir + "/scaffold-truncated/tmp\n" 
                        + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                        + "spades-res -s mapped_reads_1.fastq "
                        + "--trusted-contigs ../scaffold.fasta -k " + spadesKmerlen
                        + " --only-assembler\n"
                        + "cd spades-res\n"
                        + "samtools faidx scaffolds.fasta\n";
                }
            }
            else {
                if (spadesKmerlen.equals("default")) {
                    cmd = "cd " + outputDir + "/scaffold-truncated/tmp\n" 
                        + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                        + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq "
                        + "--trusted-contigs ../scaffold.fasta"
                        + " --only-assembler\n"
                        + "cd spades-res\n"
                        + "samtools faidx scaffolds.fasta\n";
                }
                else {
                    cmd = "cd " + outputDir + "/scaffold-truncated/tmp\n" 
                            + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                            + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq "
                            + "--trusted-contigs ../scaffold.fasta -k " + spadesKmerlen
                            + " --only-assembler\n"
                            + "cd spades-res\n"
                            + "samtools faidx scaffolds.fasta\n";
                }
            }

            FileWriter shellFileWriter = new FileWriter(outputDir + "/scaffold-truncated/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/scaffold-truncated/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log-assembly.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private static int getScaffoldFromScaffolds() {    
        int maxScaffoldLength = 0;
        String maxScaffoldId = "";
        
        File scaffoldFile = new File(outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta.fai"));
                String str = br.readLine();
                str = str.trim();
                String[] results = str.split("\t");
                maxScaffoldId = results[0].trim();
                maxScaffoldLength = Integer.parseInt(results[1].trim());            
                br.close();
                
                if (!maxScaffoldId.equals("") && maxScaffoldLength != 0)
                {
                    String cmd = "";
                    
                    cmd = "cd " + outputDir + "/scaffold-truncated/tmp/spades-res\n"
                            + "bash filterbyname.sh "
                            + "in=scaffolds.fasta "
                            + "out=scaffold.fasta names=" + maxScaffoldId
                            + " include=t\n"
                            + "samtools faidx scaffold.fasta\n";
                    
                    FileWriter shellFileWriter = new FileWriter(outputDir + "/scaffold-truncated/run.sh");
                    shellFileWriter.write("#!/bin/bash\n");
                    shellFileWriter.write(cmd);
                    shellFileWriter.close();
        
                    ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/scaffold-truncated/run.sh");
                    builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                    Process process = builder.start();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                    while (reader.readLine() != null) {
                    }
                    process.waitFor();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        //for the first time, no spades result
        else {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(outputDir + "/scaffold-truncated/scaffold.fasta.fai"));
                String str = br.readLine();
                str = str.trim();
                String[] results = str.split("\t");
                maxScaffoldLength = Integer.parseInt(results[1].trim());            
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        return maxScaffoldLength;
    }
    
    private static int getScaffoldLenFromTruncatedExtend() {
        int scaffoldLength = 0;
        File scaffoldFile = new File(outputDir + "/scaffold-truncated/scaffold.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(outputDir + "/scaffold-truncated/scaffold.fasta.fai"));
                String str = br.readLine();
                str = str.trim();
                String[] results = str.split("\t");
                scaffoldLength = Integer.parseInt(results[1].trim());            
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        else {
            //for the first time, no truncated scaffold
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
                String str = br.readLine();
                str = str.trim();
                String[] results = str.split("\t");
                scaffoldLength = Integer.parseInt(results[1].trim());            
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return scaffoldLength;
    }
    
    private static void updateCurrentScaffold() {
        File scaffoldFile = new File(outputDir + "/scaffold-truncated/scaffold.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            String cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm scaffold.fasta*\n"
                        + "cd scaffold-truncated\n" 
                        + "cp scaffold.fasta " + outputDir + "\n"
                        + "cp scaffold.fasta.fai " + outputDir + "\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }       
    }
    
    private static void createBedForTruncatedScaffold(int truncatedLen) {
        BufferedReader br = null;
        BufferedWriter bw = null;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
            bw = new BufferedWriter(new FileWriter(outputDir + "/scaffold-truncated.bed"));
            String str = br.readLine();
            String[] results = str.split("\t");
            int scaffoldLength = Integer.parseInt(results[1].trim());
            if (scaffoldLength >= truncatedLen) {
                String scaffoldId = results[0].trim();
                bw.write(scaffoldId + "\t" + truncatedLen + "\t" + (scaffoldLength - truncatedLen) + "\n");
            }
            br.close();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void growScaffoldWithAssembly() {
        int iteration = 1;
        boolean extendContig = true;        
        int prevLength = 0;
        
        while (extendContig) 
        {
            int currentLength = getScaffoldFromScaffolds();
            if (currentLength > 300000) {
                extendContig = false;
            }
            else if (currentLength > prevLength)
            {
                prevLength = currentLength;
                createBed();
                runAlignment();
                getMappedReads();
                runSpades();
                
                File scaffoldFile = new File(outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta");
                if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) {
                    iteration++;
                    extendContig = true;
                } else {
                    extendContig = false;
                }
            }
            else {
                extendContig = false;
            }
            
            if (extendContig && iteration > 2000) {
                extendContig = false;
            }
        }        
    }
    
    private static void getTruncatedScaffoldAndExtend(int currentLength) {
        String cmd = "";
        
        //truncate 300bp
        createBedForTruncatedScaffold(300);
        
        try {
            cmd = "cd " + outputDir + "\n"
                    + "rm -r scaffold-truncated\n"
                    + "mkdir scaffold-truncated\n" 
                    + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                    + "cd scaffold-truncated\n"
                    + "samtools faidx scaffold.fasta\n";

            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    
        growScaffoldWithAssembly();      
        
        int lengthFromGrowingTruncatedScaffold = 0;
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 500bp
            createBedForTruncatedScaffold(500);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 700bp
            createBedForTruncatedScaffold(700);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 1000bp
            createBedForTruncatedScaffold(1000);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 1300bp
            createBedForTruncatedScaffold(1300);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 1500bp
            createBedForTruncatedScaffold(1500);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 1700bp
            createBedForTruncatedScaffold(1700);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
        
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend();
        if (lengthFromGrowingTruncatedScaffold > currentLength) {
            return;
        }
        else 
        {
            //truncate 2000bp
            createBedForTruncatedScaffold(2000);
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            
            growScaffoldWithAssembly();
        }
    }
    
    private static boolean checkCircularity() 
    {
        String cmd = "";
        try {
            cmd = "cd " + outputDir + "\n" 
                    + "rm -r blastn-*\n";

            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
       
        //create query and subject fasta
        BufferedReader br = null;
        BufferedWriter bw1 = null;
        BufferedWriter bw2 = null;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
            bw1 = new BufferedWriter(new FileWriter(outputDir + "/blastn-subject-1stround.bed"));
            bw2 = new BufferedWriter(new FileWriter(outputDir + "/blastn-query-1stround.bed"));
            String str = br.readLine();
            String[] results = str.split("\t");
            int scaffoldLength = Integer.parseInt(results[1].trim());
            int intReadLen = (int)avgReadLen;
            if (scaffoldLength > (intReadLen*2)) {
                String scaffoldId = results[0].trim();
                bw1.write(scaffoldId + "\t" + 0 + "\t" + (scaffoldLength - (intReadLen*2)) + "\n");
                bw2.write(scaffoldId + "\t" + (scaffoldLength - (intReadLen*2)) + "\t" + (scaffoldLength-intReadLen) + "\n");
            }
            br.close();
            bw1.close();
            bw2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        cmd = "";
        try {
            cmd = "cd " + outputDir + "\n" 
                    + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-1stround.bed -fo blastn-subject-1stround.fasta\n"
                    + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-1stround.bed -fo blastn-query-1stround.fasta\n"
                    + "samtools faidx blastn-subject-1stround.fasta\n"
                    + "samtools faidx blastn-query-1stround.fasta\n"
                    + "makeblastdb -in blastn-subject-1stround.fasta -dbtype nucl\n"
                    + "blastn -query blastn-query-1stround.fasta -db blastn-subject-1stround.fasta -num_threads 16 -outfmt '7' -out blastn-res-1stround.txt\n";

            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        //parse blastn result
        BufferedWriter bwCircularOutputLog = null;
        BufferedWriter bwOutputLog = null;
        boolean isCircular = false;
        int subjectStart = 0;
        int queryStart = 0;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/blastn-res-1stround.txt"));
            bwCircularOutputLog = new BufferedWriter(new FileWriter(outputDir + "/circularity-output-log.txt", true));
            bwOutputLog = new BufferedWriter(new FileWriter(outputDir + "/output-log.txt", true));
            String str = "";
            String[] results;
            while ((str = br.readLine()) != null && !isCircular) 
            {
                if (!str.startsWith("#")) 
                {
                    str = str.trim();
                    results = str.split("\t");
                    double percentIden = Double.parseDouble(results[2].trim());
                    int alignmentLen = Integer.parseInt(results[3].trim());
                    subjectStart = Integer.parseInt(results[8].trim());
                    queryStart = Integer.parseInt(results[6].trim());
                    int minAlignmentLength = (int)(avgReadLen*0.95);

                    if (percentIden >= minIdentityCircular
                        && alignmentLen >= minAlignmentLength) {
                        bwCircularOutputLog.write("Scaffold seems circular. Scaffold position " + results[6]
                            + " to " + results[7] + " mapped to position " + results[8] + " to " + results[9] 
                            + " with " + percentIden + "% identity and " + alignmentLen + " alignment length.\n");
                        bwOutputLog.write("Scaffold seems circular. Scaffold position " + results[6]
                            + " to " + results[7] + " mapped to position " + results[8] + " to " + results[9] 
                            + " with " + percentIden + "% identity and " + alignmentLen + " alignment length.\n");
                        isCircular = true;
                    }
                    else {
                        bwCircularOutputLog.write("Scaffold does not seem to be circular. Scaffold position " + results[6]
                            + " to " + results[7] + " mapped to position " + results[8] + " to " + results[9]
                            + " with " + percentIden + "% identity and " + alignmentLen + " alignment length.\n");
                        bwOutputLog.write("Scaffold does not seem to be circular. Scaffold position " + results[6]
                            + " to " + results[7] + " mapped to position " + results[8] + " to " + results[9]
                            + " with " + percentIden + "% identity and " + alignmentLen + " alignment length.\n");
                    }
                }                
            }
            br.close();
            bwOutputLog.close();
            bwCircularOutputLog.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        //blastn subject query 2nd round
        if (isCircular) 
        {
            try {
                br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
                bw1 = new BufferedWriter(new FileWriter(outputDir + "/blastn-subject-2ndround.bed"));
                bw2 = new BufferedWriter(new FileWriter(outputDir + "/blastn-query-2ndround.bed"));
                String str = br.readLine();
                String[] results = str.split("\t");
                int scaffoldLength = Integer.parseInt(results[1].trim());
                int intReadLen = (int)avgReadLen;
                if (scaffoldLength > (intReadLen*2)) {
                    String scaffoldId = results[0].trim();
                    bw1.write(scaffoldId + "\t" + 0 + "\t" + subjectStart + "\n");
                    bw2.write(scaffoldId + "\t" + (scaffoldLength - (intReadLen*2) - subjectStart) + "\t" + (scaffoldLength-(intReadLen*2)) + "\n");
                }
                br.close();
                bw1.close();
                bw2.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            
            cmd = "";
            try {
                cmd = "cd " + outputDir + "\n" 
                        + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-2ndround.bed -fo blastn-subject-2ndround.fasta\n"
                        + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-2ndround.bed -fo blastn-query-2ndround.fasta\n"
                        + "samtools faidx blastn-subject-2ndround.fasta\n"
                        + "samtools faidx blastn-query-2ndround.fasta\n"
                        + "makeblastdb -in blastn-subject-2ndround.fasta -dbtype nucl\n"
                        + "blastn -query blastn-query-2ndround.fasta -db blastn-subject-2ndround.fasta -num_threads 16 -outfmt '7' -out blastn-res-2ndround.txt\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        return isCircular;
    }
    
    private static void runAlignmentGetCoverage() {
        String cmd = "";
        try {
            cmd = "cd " + outputDir + "\n" 
                    + "rm -r bowtie2*\n"
                    + "rm -r samtools*\n"
                    + "bowtie2-build scaffold.fasta bowtie2-index\n"
                    + "bowtie2 -x bowtie2-index "
                    + "-1 " + read1
                    + " -2 " + read2
                    + " | samtools view -bS - | samtools view -h -F 0x04 -b - | "
                    + "samtools sort - -o bowtie2-mapped.bam\n"
                    + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n";

            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private static List<Integer> percentile(List<Integer> values, List<Double> percentiles) { 
        List<Integer> percentileResults = new ArrayList<Integer>();
        Collections.sort(values); 
        for (double percentile:percentiles) {
            int index = (int) Math.ceil((percentile / 100.00) * values.size()); 
            percentileResults.add(values.get(index - 1)); 
        }
        return percentileResults;
    }
    
    private static List<Integer> readCoverageGetPercentile() {
        List<Integer> coverages = new ArrayList<Integer>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/samtools-coverage.txt"));
            
            String str = "";
            String[] results;
            
            while ((str = br.readLine()) != null) {
                str = str.trim();
                results = str.split("\t");
                int coverage = Integer.parseInt(results[2]);
                if (coverage != 0) {
                    coverages.add(coverage);
                }                
            }
            br.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        List<Double> percentiles = new ArrayList<Double>();
        percentiles.add(15.00);
        percentiles.add(85.00);
        List<Integer> percentileResults = percentile(coverages, percentiles);
        
        coverages = null;
        
        return percentileResults;
    }
    
    private static boolean writeCoverageQuantile (int quantile15Percent, int quantile85Percent) 
    {        
        boolean hasSuspiciousRegion = false;
        int scaffoldLength = 0;
        String scaffoldId = "";
        List<Integer> suspiciousStarts = new ArrayList<Integer>();
        List<Integer> suspiciousEnds = new ArrayList<Integer>();
        
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
            String str = br.readLine();
            str = str.trim();
            String[] results = str.split("\t");
            scaffoldLength = Integer.parseInt(results[1].trim()); 
            scaffoldId = results[0].trim();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        try {                       
            String str = "";
            String[] results;
            int coverage = 0;
            int startBase = 0;
            int currentBase = 0;
            boolean lowHigh = false;
            
            br = new BufferedReader(new FileReader(outputDir + "/samtools-coverage.txt"));          
            
            while ((str = br.readLine()) != null) {
                str = str.trim();
                results = str.split("\t");
                coverage = Integer.parseInt(results[2]);
                currentBase = Integer.parseInt(results[1]);
                if (coverage >= quantile15Percent && coverage <= quantile85Percent) {
                    if (lowHigh) {
                       if ((currentBase - startBase) > minSuspiciousLen) {
                           suspiciousStarts.add(startBase);
                           suspiciousEnds.add(currentBase-1);
                       }
                    }
                    lowHigh = false;
                }
                else 
                {
                    if (startBase == 0 || !lowHigh) {
                        startBase = currentBase;
                    }
                    lowHigh = true;
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }  
        BufferedWriter bwCoverageOutputLog = null;
        if (suspiciousStarts.isEmpty()) {
            hasSuspiciousRegion = false;
            try {
                bwCoverageOutputLog = new BufferedWriter(new FileWriter(outputDir + "/suspicious-regions.log", true));
                bwCoverageOutputLog.write(scaffoldId + " has no suspicious region\n");
                bwCoverageOutputLog.close();
            } catch (Exception e) {
                e.printStackTrace();
            }  
        }
        else 
        {
            hasSuspiciousRegion = true;
            int maxStartBase = 1;
            int maxLength = 0;
            int start = 0;

            try {
                bwCoverageOutputLog = new BufferedWriter(new FileWriter(outputDir + "/suspicious-regions.log", true));
                bwCoverageOutputLog.write(scaffoldId + " has " + suspiciousStarts.size() + " suspicious regions\n");
                               
                //get start and end of longest non-suspicious region  
                for (int i = 0; i < suspiciousStarts.size(); i++) {
                    start = suspiciousStarts.get(i);
                    if (maxLength == 0) {
                        maxLength = start - 1;                
                    }
                    else {
                        if ((start - suspiciousEnds.get(i-1)) > maxLength) {
                            maxLength = start - suspiciousEnds.get(i-1);
                            maxStartBase = suspiciousEnds.get(i-1) + 1;
                        }
                    }
                }
                
                if ((scaffoldLength - suspiciousEnds.get(suspiciousEnds.size()-1)) > maxLength) {
                    maxLength = scaffoldLength - suspiciousEnds.get(suspiciousEnds.size()-1);
                    maxStartBase = suspiciousEnds.get(suspiciousEnds.size()-1) + 1;
                }
                
                bwCoverageOutputLog.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            //write bed file
            BufferedWriter bw = null;
            try {
                bw = new BufferedWriter(new FileWriter(outputDir + "/longest-non-suspicious.bed"));
                bw.write(scaffoldId + "\t" + maxStartBase + "\t" + (maxStartBase + maxLength - 1) + "\n");
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            
            //get fasta from bed file
            try {
                String cmd = "cd " + outputDir + "\n"
                        + "rm -r longest-non-suspicious.fasta*\n"
                        + "bedtools getfasta -fi scaffold.fasta -bed longest-non-suspicious.bed -fo longest-non-suspicious.fasta\n"
                        + "samtools faidx longest-non-suspicious.fasta\n";
    
                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();
    
                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        
        return hasSuspiciousRegion;
    }
    
    private static void updateScaffoldWithLongest() {      
        File scaffoldFile = new File(outputDir + "/scaffold.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            String cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm scaffold.fasta*\n" 
                        + "cp longest-non-suspicious.fasta scaffold.fasta\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        
    }
    
    private static void copyScaffoldForAssembly() {
        File scaffoldFile = new File(outputDir + "/scaffold.fasta");
        if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
        {
            String cmd = "";
            try {
                cmd = "cd " + outputDir + "\n"
                        + "rm -r scaffold-truncated\n"
                        + "mkdir scaffold-truncated\n"
                        + "cp scaffold.fasta scaffold-truncated/scaffold.fasta\n"
                        + "cd scaffold-truncated\n"
                        + "samtools faidx scaffold.fasta\n";

                FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
                shellFileWriter.write("#!/bin/bash\n");
                shellFileWriter.write(cmd);
                shellFileWriter.close();

                ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
                builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
                Process process = builder.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                while (reader.readLine() != null) {
                }
                process.waitFor();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }
    
    private static boolean checkCoverage() {
        runAlignmentGetCoverage();
        List<Integer> percentiles = readCoverageGetPercentile();
        boolean hasSuspiciousRegion = writeCoverageQuantile(percentiles.get(0), percentiles.get(1));          
        
        BufferedWriter bwOutLog = null;
        try {
            bwOutLog = new BufferedWriter(new FileWriter(outputDir + "/output-log.txt", true));
            if (hasSuspiciousRegion) {                 
                updateScaffoldWithLongest();  
                bwOutLog.write("Scaffold got updated for suspicious region.\n");
                bwOutLog.close();
                copyScaffoldForAssembly();
                growScaffoldWithAssembly();                
            }
            else {
                bwOutLog.write("Scaffold didn't get updated as there is no suspicious region.\n");
                bwOutLog.close();
            }           
        } catch (Exception e) {
            e.printStackTrace();
        }
        return hasSuspiciousRegion;       
    }
    
    private static void extendOneScaffold() {
        int iteration = 1;
        String cmd = "";
        try {
            cmd = "cd " + outputDir + "\n"
                + "cp " + scaffoldFile + " scaffold.fasta\n"
                + "samtools faidx scaffold.fasta\n";
        
            FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();
    
            ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
            }
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        boolean extendContig = true;
        int prevLength = 0;
        boolean needUpdate = false;
        
        while (extendContig) {           
            int currentLength = getScaffoldLenFromTruncatedExtend();
            if (currentLength > 300000) {
                extendContig = false;
            }
            else if (currentLength > prevLength)
            {
                prevLength = currentLength;
                updateCurrentScaffold();
                if (checkCircularity())
                {
                    extendContig = false;
                }
                else
                {
                    needUpdate = checkCoverage();
                    /*if (needUpdate) {
                        //update current length
                        int newLength = getScaffoldLenFromTruncatedExtend();
                        if (newLength <= currentLength) {
                            updateCurrentScaffold();
                            prevLength = newLength;
                            getTruncatedScaffoldAndExtend(newLength);
                        }                                                
                    }*/
                    //alternate approach
                    if (needUpdate) {
                        //update current length
                        int newLength = getScaffoldLenFromTruncatedExtend();
                        updateCurrentScaffold();
                        prevLength = newLength;
                        getTruncatedScaffoldAndExtend(newLength); 
                    }
                    else {
                        getTruncatedScaffoldAndExtend(currentLength);
                    }
                    iteration++;
                }               
            }
            else {
                extendContig = false;
            }
            
            if (extendContig && iteration > 2000) {
                extendContig = false;
            }
        }
    }
    
    public static void main(String[] args) {
        parseArguments(args);
        System.out.println("Finished parsing input arguments");
        getReadLen();
        System.out.println("Started growing scaffold: " + LocalDateTime.now());      
        extendOneScaffold();
        System.out.println("Finished growing scaffold: " + LocalDateTime.now());
    }
}
