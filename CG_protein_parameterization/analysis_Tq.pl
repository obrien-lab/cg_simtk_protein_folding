#!/usr/bin/perl

my ($dt, $nsavc, $num_traj, $max_step);
my @fold = ();
print "-> Parsing info.log\n";
open(INFO, "<info.log") || die("Error: Cannot find info.log\n\n");
while(my $line = <INFO>)
{
  chomp($line);
  if($line =~ /^Number of trajectories:/)
  {
    my @str = split(/:/, $line);
    $num_traj = $str[1];
    $num_traj =~ s/^\s+|\s+$//g;
  }
  elsif($line =~ /^File save steps:/)
  {
    my @str = split(/:/, $line);
    $nsavc = $str[1];
    $nsavc =~ s/^\s+|\s+$//g;
  }
  elsif($line =~ /^Time step:/)
  {
    my @str = split(/:/, $line);
    $dt = $str[1];
    $dt =~ s/^\s+|\s+$//g;
  }
  elsif($line =~ /^\s+SIM_ID /)
  {
    my $tag = 0;
    while($line = <INFO>)
    {
      chomp($line);
      $line =~ s/^\s+|\s+$//g;
      my @str = split(/\s+/, $line);
      my $step = $str[2];
      my $if_fold = $str[4];

      if($if_fold eq 0)
      {
        print("Error: Simulation ".$str[0]." is not folded.\n");
        $tag = 1;
        push@fold, [];
      }
      else
      {
        my @str = split(/@/, $step);
        my $step = $str[1] / $nsavc;
        if($step > $max_step)
        {
          $max_step = $step;
        }
        my @array = ();
        $array[$step-1] = 1;
        push@fold, \@array;
      }
    }

    if($tag)
    {
      die();
    }
  }
}
close(INFO);
print "   Done. time_step = $dt, num_traj = $num_traj, nsteps_save = $nsavc, max_step = $max_step\n";

print "-> Calculating folding status\n";
my $time_step_ns = $dt * $nsavc / 1000;

if($num_traj ne @fold)
{
  die("Error: number of trajectories mismatched\n");
}

for(my $i = 0; $i < $num_traj; $i++)
{
  my @array = @{$fold[$i]};
  for(my $j = @array+1; $j <= $max_step; $j++)
  {
    push@array, 1;
  }
  $fold[$i] = \@array;
}

open(LOG, ">Su_vs_time.dat") || die("Error: Cannot create Su_vs_time.dat\n\n");
printf LOG ("%-15s %-10s\n", "Time(ns)", "Su");
for(my $i = 0; $i < $max_step; $i++)
{
  my $Nn = 0;
  for(my $x = 0; $x < $num_traj; $x++)
  {
    if($fold[$x]->[$i])
    {
      $Nn++;
    }
  }
  my $Su = 1 - $Nn / $num_traj;
  $Su = sprintf("%.6f", $Su);
  my $time = ($i+1) * $time_step_ns;
  $time = sprintf("%.3f", $time);
  printf LOG ("%-15s %-10s\n", $time, $Su);
}
print "   Done\n";

print "-> Fitting survival probability vs time\n";
open(M, ">fit_k.m") || die("Error: Cannot create fit_k.m\n\n");
print M "function fit_k()
clear all
fileID = fopen('Su_vs_time.dat');
C = textscan(fileID,'%f %f','HeaderLines',1);
fclose(fileID);
t=C{1};
Su=C{2};

index_1 = 1;
index_0 = find(Su == 0);
t=t(index_1(end):index_0(1));
Su=Su(index_1(end):index_0(1));

x0=[0 0];
[x,resnorm] = lsqcurvefit(\@exp_fun,x0,t,Su);

Su_fit=exp_fun(x,t);
[R,P]=corrcoef(Su, Su_fit);
disp(['R2 = ', num2str(R(1,2))]);
disp(['p = ', num2str(P(1,2))]);
k=x(1);
t0=x(2);
disp(['k_F = ' sprintf('%.3f',k)]);
disp(['t0 = ' sprintf('%.3f',t0)]);
disp(['tao_F = ' sprintf('%.2f',1/k+t0)]);

t1=t0:0.001:t(end);
Su1=exp_fun(x,t1);
figure(1)
hold on
set(gcf,'Units','centimeters','Position',[8 1.5 10 8],'paperpositionmode','auto');
plot(t,Su,'ok','LineWidth',0.1,'MarkerSize',3)
plot(t1,Su1,'-r','LineWidth',3)
set(gca, 'fontsize',12,'fontweight','normal','LineWidth',1.0,'fontname','Times')
box on
grid on
xlabel('\$t\\ (\\rm ns)\$','fontsize',14,'fontweight','b','color','k','Interpreter','latex')
ylabel('\$S_{\\rm U}\$','fontsize',14,'fontweight','b','color','k','Interpreter','latex')
legend('String',{'\$S_{\\rm U}^{\\rm clac}\$', ['\$S_{\\rm U}=e^{-\\frac{1}{' sprintf('%.2f',1/k) '}(t-'...
    sprintf('%.3f',t0) ')}\$']},...
    'Location','best','fontsize',12, 'fontweight','normal','box','off','Interpreter','latex');
saveas(gcf, 'fit_k.fig');
quit
end

function f = exp_fun(x,xdata)
f = exp(-x(1)*(xdata-x(2)));
end\n";
close(M);

my $screen = `matlab -nodisplay -nosplash -r fit_k 2>&1`;
if($screen =~ /tao_F = /)
{
  chomp($screen);
  my @str = split(/\n/, $screen);
  foreach my $line (@str)
  {
    if($line =~ /R2 =/)
    {
      $line =~ s/^\s+|\s+$//g;
      print "   $line\n";
    }
    elsif($line =~ /p =/)
    {
      $line =~ s/^\s+|\s+$//g;
      print "   $line\n";
    }
    elsif($line =~ /k_F =/)
    {
      $line =~ s/^\s+|\s+$//g;
      print "   $line\n";
    }
    elsif($line =~ /t0 =/)
    {
      $line =~ s/^\s+|\s+$//g;
      print "   $line\n";
    }
    elsif($line =~ /tao_F =/)
    {
      $line =~ s/^\s+|\s+$//g;
      print "   $line\n";
    }
  }
}
else
{
  die("Error: Failed to fit\n$screen\n");
}
print "   Done\n";
