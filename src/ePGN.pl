#!/usr/bin/perl -w

#
# ePGN.pl - Programa para estimacao de uma Probabilistic Genetic Network (PGN).
#
# Argumentos:
# 
# -i=arquivo.txt : arquivo que contem uma serie temporal de uma PGN
# -c : utilizar CoD: Coeficiente de Determinacao (default: Informacao Mutua Media -- IMM)
# -p=n : estimar considerando n preditores para o gene-alvo (default:2)
# -n=m : estimar IMM considerando o parametro n com valor m (default:0)
#
# Devolve, na saida padrao, uma tabela texto, contendo a lista dos melhores preditores para cada gene,
# ordenados em ordem descrescente de pontuacao.
# 
# Copyleft M.S.Reis, 14 de agosto 2009
#

use strict;

# Definicoes
#
my $MAX_NUM_ITER = 1000;
my $MIN_NUM_ITER = 1;
my $MAX_NOISE = 100;
my $MIN_NOISE = 0;
my $MAX_VALUE = 2;   # discretizado em 3 estados: 0, 1 ou 2
my $MIN_VALUE = 0;   # por hora somente valores nao-negativos

# Variaveis
#
my $estimador = "e";
my $num_preditores = 2;
my $nome_arq = "";
my $T = 0;    # instantes de tempo da serie temporal
my $N;        # numero de genes
my $n = 0;    # parametro utilizado pelo IMM, para cortar valores com poucas amostras

# Le os argumentos
#
foreach(@ARGV){
  if ($_ =~ /\-c/){
    $estimador = "c";
  }
  elsif ($_ =~ /\-i\=(.*)/){
    $nome_arq = $1;
  }
  elsif ($_ =~ /\-n\=(\d+)/){
    $n = $1;
    if ($1 < 0){
      print "Erro! Parametro n deve ser inteiro nao-negativo!\n";
      erroParametro();
    }
  }
  elsif ($_ =~ /\-p\=(\d+)/){
    $num_preditores = $1;
    if ($1 <= 0){
      print "Erro! Numero de preditores deve ser inteiro positivo!\n";
      erroParametro();
    }
  }
  else{
    erroParametro();
  }
}

# Abre o arquivo de serie temporal da rede (PGN)
#
if (($nome_arq eq "") or (!open(ARQ,$nome_arq))){
  print "\nNome de arquivo ($nome_arq) invalido ou arquivo inexistente!\n";
  erroParametro();
}

# Carrega o arquivo de entrada
#
my %S;    # estado dos genes nos T instantes de tempo
my @gene; # indexa os nomes dos genes em uma lista

while(<ARQ>){
  chomp $_;
  my @lista = split (" ", $_);
  $T = $#lista;
  my $chave = shift @lista;
  push @gene, $chave;
  $S{$chave} = [ @lista ];
}
close(ARQ);

$N = scalar keys %S;

($N >= $num_preditores) or die "Erro: num. de preditores ($num_preditores) maior que num. de genes ($N)!\n";

my %H;
my %E;
my $COD_Y;

foreach my $alvo (sort keys %S){

  if($estimador eq "c"){
    #
    # Estimando utilizando o CoD - MSE como medida de erro.
    # Erro do classificador verificado utilizando ressubstituicao (treina o classificador com toda a serie temporal
    # e utiliza a propria serie temporal para avaliar o erro do classificador)

    my $Y = 0;
    $COD_Y = 0;

    # Calcula a variancia (MSE da media):
    #
    foreach(my $i=1; $i <= $#{$S{$alvo}}; $i++){
      $Y += ${$S{$alvo}}[$i];
    }
    $Y = int( $Y / ($T-1) );
    foreach(my $i=1; $i <= $#{$S{$alvo}}; $i++){
      $COD_Y += ( $Y - ${$S{$alvo}}[$i] ) ** 2;
    }

    todasAsCombinacoes($N, $num_preditores, $alvo, "COD");   

    printf "\nPredicoes de $alvo (utilizando COD):\n";
    foreach my $pred (reverse sort {$E{$alvo}->{$a} <=> $E{$alvo}->{$b}} keys %{$E{$alvo}}){
      printf "COD\_$alvo\_($pred) = %3.3f\n", $E{$alvo}->{$pred};
    }
	
  }
  
  else{
    #
    # Estimando utilizando a informacao mutua media (IMM).
    #
    
    # Y = variavel aleatoria valor do gene-alvo no instante de tempo t
    # X = v.a. dos genes preditores de Y, no instante t-1
    # Pr(Y = y) = probabilidade de Y ter valor y    Pr:{0,1,2}->[0,1]
    # Entropia de Y: H(Y) = - SUM (Pr(Y = y) * log Pr(Y=y))
    # Informacao Mutua Media: E[I(X,Y)] = H(Y) - E[H(Y|X)]
    # E[H(Y|X)] = - SUM { Pr(X) * SUM[ Pr(Y|X) * log(Pr(Y|X)) ] }
    #
    my %Pr;

    # IMPORTANTE: talvez nao possa ser considerado o t=0 para Pr(Y), pois o Pr(Y|X) comeca com t=1!!!
    # (atualmente eh considerado a partir de t=1 para calcular Pr(Y))
    foreach(my $i=1; $i <= $#{$S{$alvo}}; $i++){
      !defined($Pr{$alvo}->[${$S{$alvo}}[$i]]) and $Pr{$alvo}->[${$S{$alvo}}[$i]] = 0  or $Pr{$alvo}->[${$S{$alvo}}[$i]]++;
    }
  
    $H{$alvo} = 0;

    for(my $i = $MIN_VALUE; $i <= $MAX_VALUE; $i++){
      if(!defined $Pr{$alvo}->[$i]){
	$Pr{$alvo}->[$i] = 0; 
      }
      $Pr{$alvo}->[$i] = $Pr{$alvo}->[$i] / ($T-1);
      if($Pr{$alvo}->[$i] > 0){
	$H{$alvo} = $H{$alvo} - (log($Pr{$alvo}->[$i]) * $Pr{$alvo}->[$i]);
      }
      # printf "gene: $alvo, T = $T, Pr($i) = %3.3f, \n", $Pr{$alvo}->[$i];
    }
    # printf "H = %3.3f\n\n", $H{$alvo};
    
    todasAsCombinacoes($N, $num_preditores, $alvo, "IMM");   

    printf "\nPredicoes de $alvo (utilizando IMM):\n";
    foreach my $pred (reverse sort {$E{$alvo}->{$a} cmp $E{$alvo}->{$b}} keys %{$E{$alvo}}){
      printf "IMM\_$alvo\_($pred) = %3.3f\n", $E{$alvo}->{$pred};
    }
  }  # fim do else (IMM)

}
  

# Fim do "main"
#
exit 0;

#            #
# Subrotinas #
#            #

# Subrotina que exibe a sintaxe do programa e retorna ao sistema.
#
sub erroParametro{

  die "\n\nSintaxe: ./ePGN.pl [-c | -p=n | -n=m] -i=serie_temporal.txt\n\n";

}


# Subrotina que implementa o algoritmo de Knuth para gerar todas as (s,t)-combinacoes
# Algoritmo (Algorithm L) disponivel em TACP 4.2.3.
# Recebe o gene-alvo, N e t. Supoe N >= t >= 0. 
#
sub todasAsCombinacoes{
  my ($N, $t, $alvo, $metodo) = @_;
  my $j;
  my @c;
  # L1. inicializacao
  for($j = 1; $j <= $t; $j++){
    $c[$j] = $j - 1;
  }
  $c[$t+1] = $N;
  $c[$t+2] = 0;

  while(1){
    
    # L2. "visita" a combinacao
    #
    # print "$c[$_] " foreach 1..$t;
    # print "\n";
    if($metodo eq "IMM"){
      calculaProbabilidadeIMM($alvo, @c);
    }
    else{
      calculaProbabilidadeCOD($alvo, @c);
    }

    # L3. "encontra" j
    #
    $j = 1;
    while($c[$j]+1 == $c[$j+1]){
      $c[$j] = $j-1;
      $j++;
    }
    # L4. done?
    # 
    if ($j > $t){
      return;
    }
    # L5. incrementa c_j
    #
    else{
      $c[$j]++;
    }
  }
}


# 
#
#
sub calculaProbabilidadeCOD{
  my ($alvo, @preditores) = @_;

  
  pop @preditores; # dois pops pra remover os elementos que nao sao preditores
  pop @preditores; #
  
  my %classificador;
  my $Err_Classificador = 0;


  # Treina o classificador com a serie temporal
  #
  for (my $i = 1; $i <= $T; $i++){
    my $pred = "";

    foreach (1..$#preditores){
      $pred = $pred . " " . $S{$gene[$preditores[$_]]}->[$i-1];
    }
    if(!defined($classificador{$pred}->[$S{$alvo}->[$i]])){
      $classificador{$pred}->[$S{$alvo}->[$i]] = 0;
    }
    else{
      $classificador{$pred}->[$S{$alvo}->[$i]]++;	    
    }    
  }
  
  # Verifica o erro do classificador (MSE) utilizando as mesmas amostras (ressubstituicao)
  #
  for (my $i = 1; $i <= $T; $i++){
    my $pred = "";

    foreach (1..$#preditores){
      $pred = $pred . " " . $S{$gene[$preditores[$_]]}->[$i-1];
    }
    my $mais_frequente = 0;
    my $X = $MIN_VALUE;
    for(my $j = $MIN_VALUE; $j <= $MAX_VALUE; $j++){
      if(!defined($classificador{$pred}->[$j])){
	$classificador{$pred}->[$j] = 0;
      }
      if($classificador{$pred}->[$j] >= $mais_frequente){
	$mais_frequente = $classificador{$pred}->[$j];
	$X = $j;
      }
      $Err_Classificador += ($S{$alvo}->[$i] - $X) ** 2;
    }
  }

  # computa o CoD
  #
  my $pred = "";
  foreach (1..$#preditores){
    $pred = $pred . " " . $gene[$preditores[$_]];
  }
  $E{$alvo}->{$pred} = ($COD_Y - $Err_Classificador) / $COD_Y;

}
  

#
#
#
sub calculaProbabilidadeIMM{
  my ($alvo, @preditores) = @_;
  
  pop @preditores; # dois pops pra remover os elementos que nao sao preditores
  pop @preditores; #
  
  my @estado;

  my $pred = "";
  my $sum = 0;

  # inicializa a combinacao atual (x1,x2,...,xk) com o valor minimo (X=(0,0,....,0))
  #
  foreach (1..$#preditores){
    # print "$gene[$preditores[$_]] ";
    $pred = $pred . " " . $gene[$preditores[$_]];
    $estado[$preditores[$_]] = $MIN_VALUE;
  }

  # Pr(X=(a,b)) * SUM_y (Pr(Y=y | X=(a,b)) * log (Pr(Y=y | X=(a,b))))
  #

  # 
  do{
    my $i = 1;
    while(($estado[$preditores[$i]] > $MAX_VALUE)and($i < $#preditores)){
      $estado[$preditores[$i]] = $MIN_VALUE;
      $i++;
      $estado[$preditores[$i]]++;  #  "vai um"
    }
    
    if($i < $#preditores){
    
      my @freq_alvo;
      my @freq_alvo_e_pred;
      my $freq_pred = 0;
      my $sum_y = 0;

      for(my $j = $MIN_VALUE; $j <= $MAX_VALUE; $j++){
	$freq_alvo[$j] = $freq_alvo_e_pred[$j] = 0;
      }
      
      for(my $k = 1; $k < $T; $k++){ 
	
	$freq_alvo[$S{$alvo}->[$k]]++;    # conta o Pr(Y=y)

	my $verifica_estado_pred = 1;
	
	foreach (1..$#preditores){
	  
	  if ($S{$gene[$preditores[$_]]}->[$k-1] != $estado[$preditores[$_]]){
	    $verifica_estado_pred = 0;
	  }
	}
	if($verifica_estado_pred){
	  $freq_pred++;
	  $freq_alvo_e_pred[$S{$alvo}->[$k]]++; # conta o Pr(Y=y e X=(a,b))
	}
      }
      
      # para depuracao
      #
      # print "\n(";
      # foreach(1..$#preditores){
	# printf "%d ", $estado[$preditores[$_]];
      # }
      # print ")\n";
      for(my $j = $MIN_VALUE; $j <= $MAX_VALUE; $j++){
	# printf "valor_$alvo = %d, freq_alvo = %d, freq_pred = %d, freq_alvo_e_pred = %d\n",
        #           $j, $freq_alvo[$j], $freq_pred, $freq_alvo_e_pred[$j];
	my $prob = 0;
	if($freq_pred > 0){
	  $prob = $freq_alvo_e_pred[$j] / $freq_pred;
	}
	if($prob > 0){
	  $sum_y += $prob * log($prob);
	}
	# printf "Pr(Y=$j | X=....) = %3.3f\n", $prob;
      }

      $sum -= $sum_y * ($freq_pred / ($T-1));   #  Pr(X=(a,b)) *  SUM_y (Pr(Y=y|X=(a,b) * log(....

      $estado[$preditores[1]]++;      
    }

  }while($estado[$preditores[$#preditores]] <= $MAX_VALUE);
  
  #  printf "\nPredicao de $alvo utilizando ($pred )";
  #  printf "\nE[H(Y|X)] = %3.3f\n", $sum;
  #  printf "E[I(X,Y)] = %3.3f\n", $H{$alvo} - $sum;
  #  print "\n";

  $E{$alvo}->{$pred} = $H{$alvo} - $sum;

}

