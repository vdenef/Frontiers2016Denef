#/usr/bin/ruby -w

unless ARGV.length == 7
  puts "This script modifies tab delimited file by replacing names in one column with new ones based on input file. Input arguments are (1) tab-delimited file with names to swap and replacement names (2) column number, starting with 0, containing term to search for, (3) column number containg term that will replace searched term (4) tab-delimited file with data fields to be converted (5) column matching name to search for (6) column in data file to replace (7) a file containing the pattern to match based on which to select OTUs to do the replacement for"
  exit(1)
end

puts "Setting up replacement table to replace term in column #{ARGV[1]} with #{ARGV[2]} of input file #{ARGV[0]}"

list = Hash.new
File.open(ARGV[0]) do |file|
  file.each do |line|
    line.chomp!
    list[line.split("\t")[ARGV[1].to_i]] = line.split("\t")[ARGV[2].to_i]
  end
end


match = String.new
File.open(ARGV[6]) do |file|
  file.each do |line|
    match = line.chomp!
  end
end

puts "replacing taxonomy for OTUs matching #{match} and writing to outfile"

outfile1 = File.new("#{ARGV[3]}"'.renamed',"w")
File.open(ARGV[3]) do |file|
  file.each do |line|
    line.chomp!
    i = ARGV[4].to_i
    j = ARGV[5].to_i	
    temp = line.split("\t")
    if (temp[j] =~ /#{match}/)
        temp[j] = list[temp[i]]
	puts temp
    end
    outfile1.puts "#{temp.join("\t")}"
  end
end

outfile1.close

