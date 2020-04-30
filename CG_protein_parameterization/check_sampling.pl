#!/usr/bin/perl
use POSIX;

my $log_file = $ARGV[0];
my $skip_step = $ARGV[1];
my $type = $ARGV[2];

if($skip_step eq "")
{
  $skip_step = 0;
}

if($type eq "")
{
  $type = 0;
}

my @file_list = split(/\s+/, $log_file);

my $step = 0;
my $n_windows;
my @counter_all = ();

open(DAT, "<".$file_list[0]) || die("Error: Cannot find ".$file_list[0]."\n\n");
while(my $line = <DAT>)
{
  if($line =~ /^Number of windows:/)
  {
    chomp($line);
    my @str = split(/:/, $line);
    $n_windows = $str[1];
    $n_windows =~ s/^\s+|\s+$//g;
    for(my $i = 0; $i < $n_windows; $i++)
    {
      push@counter_all, 0;
    }
  }
  elsif($line =~ /^Temperatures: /)
  {
    chomp($line);
    my @str = split(/\s+/, $line);
    @temps = @str[1...$#str];
  }
}
close(DAT);

my @end_points = ();
my @idx_map_list = ();

foreach my $file (@file_list)
{
  if($file =~ /_r_/)
  {
    my @str = split(/_r_/, $file);
    @str = split(/.log/, $str[$#str]);
    push@end_points, $str[0];
  }
}
push@end_points, 999999999999999;

mkdir("analysis");
mkdir("analysis/sampling");
open(OUT, ">analysis/sampling/replica.log");

my @tag = ();
my @counter = ();
for(my $i = 0; $i < $n_windows; $i++)
{
  push@counter, 0;
  push@tag, -1;
}

for(my $k = 0; $k < @file_list; $k++)
{
  my $file = $file_list[$k];
  my $end_point = $end_points[$k];
  my @idx_map = ();
  for(my $i = 1; $i <= $n_windows; $i++)
  {
    $idx_map[$i] = $i;
  }
  for(my $kk = $k - 1; $kk >= 0; $kk--)
  {
    my @idx_map_0 = @{$idx_map_list[$kk]};
    for(my $i = 1; $i <= $n_windows; $i++)
    {
      my $ii = $idx_map_0[$idx_map[$i]];
      $idx_map[$i] = $ii;
    }
  }
  #print "@idx_map\n";
  
  open(DAT, "<".$file) || die("Error: Cannot find ".$file."\n\n");
  my $end_line = "";
  while(my $line = <DAT>)
  {
    if($line =~ /^PROD [0-9]+:/)
    {
      chomp($line);
      my @str = split(/:/, $line);
      @str = split(/\s+/, $str[0]);
      $step = $str[1];
      $step =~ s/^\s+|\s+$//g;
      if($step > $skip_step && $step < $end_point)
      {
        for(my $i = 1; $i <= $n_windows; $i++)
        {
          if($line =~ /: \Q$i\E\|/)
          {
            if($tag[$idx_map[$i]-1] eq -1 && $counter[$idx_map[$i]-1] eq 0)
            {
              $tag[$idx_map[$i]-1] = 0;
            }
            if($tag[$idx_map[$i]-1] eq 0 && $counter[$idx_map[$i]-1] % 2 ne 0)
            {
              $counter[$idx_map[$i]-1]++;
            }
            elsif($tag[$idx_map[$i]-1] eq 1 && $counter[$idx_map[$i]-1] % 2 eq 0)
            {
              $counter[$idx_map[$i]-1]++;
            }
          }
          elsif($line =~ / \Q$i\E\| /)
          {
            if($tag[$idx_map[$i]-1] eq -1 && $counter[$idx_map[$i]-1] eq 0)
            {
              $tag[$idx_map[$i]-1] = 1;
            }
            if($tag[$idx_map[$i]-1] eq 0 && $counter[$idx_map[$i]-1] % 2 eq 0)
            {
              $counter[$idx_map[$i]-1]++;
            }
            elsif($tag[$idx_map[$i]-1] eq 1 && $counter[$idx_map[$i]-1] % 2 ne 0)
            {
              $counter[$idx_map[$i]-1]++;
            }
          }
        }
        my @str = split(/: /, $line);
        @str = split(/\|+ /, $str[1]);
        print OUT "$step ";
        for(my $i = 1; $i <= $n_windows; $i++)
        {
          print OUT $idx_map[$str[$i-1]]." ";
        }
        print OUT "\n";
        
        $end_line = $line;
      }
    }
  }
  
  my @str = split(/: /, $end_line);
  @str = split(/\|+ /, $str[1]);
  my @idx_map_1 = ();
  for(my $i = 1; $i <= $n_windows; $i++)
  {
    $idx_map_1[$i] = $str[$i-1];
  }
  push@idx_map_list, \@idx_map_1;
  close(DAT);
}
close(OUT);

for(my $i = 0; $i < $n_windows; $i++)
{
  $counter_all[$i] = int($counter[$i]/2);
}

for(my $t = 1; $t <= @temps; $t++)
{
  printf("%6s ", "# $t");
}
print("\n");

foreach my $c (@counter_all)
{
  printf("%6s ", $c);
}
print("\n");

chdir("analysis/sampling");

Ep_distribution_plot(\@temps, $skip_step);

if($type eq 1)
{
  high_low_exch_plot(\@temps, \@counter_all);
  Cv_vs_step_plot(\@temps, $skip_step, $step, 10000);
}
elsif($type eq 2)
{
  high_low_exch_plot(\@temps, \@counter_all);
}
elsif($type eq 3)
{
  Cv_vs_step_plot(\@temps, $skip_step, $step, 10000);
}
######################################################
sub high_low_exch_plot
{
  my @temps = @{$_[0]};
  my @counter_all = @{$_[1]};
  my $nsim_temps = @temps;
  open(MAT, ">plot_exchange.m");
  print MAT "data = load('replica.log');
n_window = $nsim_temps;
T_list = [";
  foreach my $t (@temps)
  {
  	print MAT $t." ";
  }
  print MAT "];
counter = [";
  foreach my $c (@counter_all)
  {
  	print MAT $c." ";
  }
  print MAT "];

for i = 1:size(data, 1)
    for j = 1:n_window
        replica_idx = data(i, j+1);
        D(i,replica_idx) = j;
    end
end

figure(1)
set(gcf,'Units','centimeters','Position',[8 1.5 35 25],'paperpositionmode','auto');
x = data(1,1):data(end,1);
for j=1:n_window
    subplot(3,6,j)
    plot(x, D(:,j), '-k')
    set(gca, 'YTick', 1:length(T_list), 'YTickLabel', T_list, 'YLim', [1 length(T_list)], ...
        'YGrid', 'on')
    title(['# ' num2str(j) ': ' num2str(counter(j))])
end

saveas(gcf, 'exchange.svg')
quit";
  close(MAT);
  print("-> Going to plot high to low exchange...\n");
  system("matlab -nodisplay -r plot_exchange > /dev/null");
  print("   Done\n");
}

#######################################################
sub Cv_vs_step_plot
{
  my @temps = @{$_[0]};
  my $skip_step = $_[1];
  my $tot_step = $_[2];
  my $ds = $_[3];
  my $max_temp = $temps[$#temps];
  $max_temp += 20;
  my $min_temp = $temps[0];
  $min_temp -= 20;
  my $nsim_temps = @temps;

  print("-> Going to calculate Cv vs different steps...\n");
  for(my $i = 1; $i <= POSIX::ceil(($tot_step-$skip_step)/$ds); $i++)
  {
  	open(IN, ">myreaders.py");
  	print IN "def ReadLastColumn(filename):
    QList = []
    results = []
    n = 1
    f = open(filename, \"r\")
    for line in f:
        if n > ".($skip_step+$ds*$i)." or n < ".($skip_step+$ds*($i-1)+1).":
            n = n + 1
            continue
        QList.append(float(line.strip()))
        n += 1
    f.close()
    results.append(QList)
    return results";
    close(IN);
    open(CNTRL, ">wham_cv.xml");
    print CNTRL "<?xml version=\"1.0\" ?>
<WhamSpec>
    <General>
        <Coordinates>
            <Coordinate name=\"V\" />
        </Coordinates>
        <CoordinateFileReaders pythonModule=\"myreaders\" returnsTime=\"false\">
            <Reader name=\"ReadLastColumn\">
                <ReturnList name=\"V\" />
            </Reader>
        </CoordinateFileReaders>
        <Binnings>
            <Binning name=\"V\">
                <Interval>1</Interval>
            </Binning>
        </Binnings>
        <Parameters>
            <Parameter name=\"kB\">0.0019872041</Parameter>
        </Parameters>
    </General>
    <Trajectories>\n";
    for(my $j=1;$j<=$nsim_temps;$j++)
    {
      print CNTRL " " x 8 . "<Trajectory T=\"" . $temps[$j-1] . "\">
            <EnergyFunction>V</EnergyFunction>
            <CoordinateFiles>
                <CoordinateFile>../../aa${j}/ene_1.log</CoordinateFile>
            </CoordinateFiles>
        </Trajectory>\n";
    }
    print CNTRL "    </Trajectories>
    <Jobs>
        <HeatCapacity outFile=\"cv/cv_$i.dat\">
            <EnergyFunction>V</EnergyFunction>
            <Temperatures>${min_temp}:0.1:${max_temp}</Temperatures>
        </HeatCapacity>
    </Jobs>
</WhamSpec>";
    close(CNTRL);
    system("pywham.py wham_cv.xml > /dev/null 2>&1");
  }

  open(MAT, ">plot_CV_sampling.m");
  print MAT "D = {};
Tm = [];
T_max = -9999999;
T_min = 9999999;
legend_str = {};
n=1.1;
for i=1:".POSIX::ceil(($tot_step-$skip_step)/$ds)."
    C = load(['cv/cv_' num2str(i) '.dat']);
    D{i} = C;
    [m,j]=max(C(:,2));
    Tm(i,:)=[C(j,1), m];
    T_max = max([T_max; C(:,1)]);
    T_min = min([T_min; C(:,1)]);
    legend_str{i} = ['\$', num2str(".$skip_step."+(i-1)*".$ds."+1) , '-', num2str(".$skip_step."+i*".$ds."), '\$'];
end
T_min = floor(T_min);
T_max = ceil(T_max);

figure(1)
hold on
set(gcf,'Units','centimeters','Position',[8 1.5 10 8],'paperpositionmode','auto');
for i=1:".POSIX::ceil(($tot_step-$skip_step)/$ds)."
    plot(D{i}(:,1),D{i}(:,2),'-','LineWidth',1.5)
end
set(gca, 'fontsize',10,'fontweight','normal','LineWidth',1.0,'fontname','Nimbus Roman No9 L')
axis([T_min T_max 0 max(Tm(:,2))+1])
box on
grid on
xlabel('\$T\\ (\\rm K)\$','fontsize',12,'color','k','Interpreter','latex')
ylabel('\$C_{\\rm {V}}\\ (\\rm {kcal/mol/K})\$','fontsize',12,'color','k','Interpreter','latex')
h = legend('String',legend_str,'Location','best','fontsize',10,'box','off','Interpreter','latex');

for i=1:".POSIX::ceil(($tot_step-$skip_step)/$ds)."
    plot([Tm(i,1) Tm(i,1)],[0 Tm(i,2)],'--k','LineWidth',1.3)
end
saveas(gcf, 'cv_sampling.svg');
quit";
  close(MAT);
  system("matlab -nodisplay -r plot_CV_sampling > /dev/null");
  print("   Done\n");
  
  open(IN, ">myreaders.py");
  print IN "def ReadLastColumn(filename):
  QList = []
  results = []
  n = 1
  f = open(filename, \"r\")
  for line in f:
      if n < ".($skip_step+1).":
          n = n + 1
          continue
      QList.append(float(line.strip()))
      n += 1
  f.close()
  results.append(QList)
  return results";
  close(IN);
  open(CNTRL, ">wham_cv.xml");
  print CNTRL "<?xml version=\"1.0\" ?>
<WhamSpec>
    <General>
        <Coordinates>
            <Coordinate name=\"V\" />
        </Coordinates>
        <CoordinateFileReaders pythonModule=\"myreaders\" returnsTime=\"false\">
            <Reader name=\"ReadLastColumn\">
                <ReturnList name=\"V\" />
            </Reader>
        </CoordinateFileReaders>
        <Binnings>
            <Binning name=\"V\">
                <Interval>1</Interval>
            </Binning>
        </Binnings>
        <Parameters>
            <Parameter name=\"kB\">0.0019872041</Parameter>
        </Parameters>
    </General>
    <Trajectories>\n";
  for(my $j=1;$j<=$nsim_temps;$j++)
  {
    print CNTRL " " x 8 . "<Trajectory T=\"" . $temps[$j-1] . "\">
            <EnergyFunction>V</EnergyFunction>
            <CoordinateFiles>
                <CoordinateFile>../../aa${j}/ene_1.log</CoordinateFile>
            </CoordinateFiles>
        </Trajectory>\n";
  }
  print CNTRL "    </Trajectories>
    <Jobs>
        <HeatCapacity outFile=\"cv/cv_all.dat\">
            <EnergyFunction>V</EnergyFunction>
            <Temperatures>${min_temp}:0.1:${max_temp}</Temperatures>
        </HeatCapacity>
    </Jobs>
</WhamSpec>";
    close(CNTRL);
    system("pywham.py wham_cv.xml > /dev/null 2>&1");

  print("-> ploting CV with windows...\n");
  open(MAT, ">plot_CV_windows.m");
  print MAT "D = {};
Tm = [];
T_max = -9999999;
T_min = 9999999;
legend_str = {};
dir='cv/cv_all.dat';
windows=[@temps];
i=1;
C = load(dir);
D{i} = C;
[m,j]=max(C(:,2));
Tm(i,:)=[C(j,1), m];
T_max = max([T_max; C(:,1)]);
T_min = min([T_min; C(:,1)]);
T_min = floor(T_min);
T_max = ceil(T_max);

figure(1)
hold on
set(gcf,'Units','centimeters','Position',[8 1.5 10 8],'paperpositionmode','auto');
for i=1
    plot(D{i}(:,1),D{i}(:,2),'-','LineWidth',1.5)
end
set(gca, 'fontsize',10,'fontweight','normal','LineWidth',1.0,'fontname','Nimbus Roman No9 L')
axis([T_min T_max 0 max(Tm(:,2))+1])
box on
grid on
xlabel('\$T\\ (\\rm K)\$','fontsize',12,'color','k','Interpreter','latex')
ylabel('\$C_{\\rm {V}}\\ (\\rm {kcal/mol/K})\$','fontsize',12,'color','k','Interpreter','latex')

plot([Tm(i,1) Tm(i,1)],[0 Tm(i,2)],'--k','LineWidth',1.3)

for i=1:length(windows)
    idx = find(D{1}(:,1) == windows(i));
    plot(windows(i), D{1}(idx,2), 'ok')
end

saveas(gcf, 'cv_windows.svg');
quit";
  close(MAT);
  system("matlab -nodisplay -r plot_CV_windows > /dev/null");
  print("   Done\n");
}

#######################################################
sub Ep_distribution_plot
{
  my @temps = @{$_[0]};
  my $skip_step = $_[1];
  
  print("-> ploting Ep distribution...\n");
  open(MAT, ">plot_Ep_distribution.m");
  print MAT "clear all
n_skip = $skip_step;
T_windows = [@temps];
legend_str = {};
for i=1:length(T_windows)
    legend_str{i} = ['T = ' num2str(T_windows(i)) ' K'];
end

figure(1)
set(gcf,'Units','centimeters','Position',[8 1.5 12 10],'paperpositionmode','auto');
C = colormap(jet(length(T_windows)));
hold on
for i=1:length(T_windows)
    D = [];
    D = load(['../../aa' num2str(i) '/ene_1.log']);
    [N, edges] = histcounts(D(n_skip+1:end),50,'Normalization','probability');
    x = (edges(1:end-1)+edges(2:end))/2;
    patch(x,N,C(i,:),'FaceAlpha',0.4,'EdgeColor','k','LineWidth',1.0)
end
set(gca, 'fontsize',11,'fontweight','normal','LineWidth',1.5,'fontname','Nimbus Roman No9 L')
box on
grid on
xlabel('Potential Energy (kcal/mol)','fontsize',13,'color','k','Interpreter','latex')
ylabel('Probability','fontsize',13,'color','k','Interpreter','latex')
h = legend('String',legend_str,'Location','best','fontsize',6,'box','off','Interpreter','latex');
saveas(gcf, 'Ep_distribution.svg');
quit";
  close(MAT);
  system("matlab -nodisplay -r plot_Ep_distribution > /dev/null");
  print("   Done\n");
}
