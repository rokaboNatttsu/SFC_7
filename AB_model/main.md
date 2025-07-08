- [TODO](#todo)
  - [バグ](#バグ)
- [1. エージェント番号](#1-エージェント番号)
- [2. モデルの式とアルゴリズム](#2-モデルの式とアルゴリズム)
- [拡張の方向性](#拡張の方向性)
  - [アルゴリズム](#アルゴリズム)
  - [行動の仮定](#行動の仮定)

# TODO

jlファイルの変更をこちらに反映させる作業をまだ行っていない

## バグ
- 倒産をきっかけに不具合が起こっている
  - 株式発行量、預金残高、借入金残高は、ストックもフローも、経済主体感の整合性は取れていると思われる
- 取引の整合性が取れていない
- 資金フローの整合性が取れていない
- 

# 1. エージェント番号
|      | インデックス | 総数 |
| :--- | :----------: | :--: |
| 家計 |     $j$      | $J$  |
| 企業 |     $o$      | $O$  |
| 銀行 |     $n$      | $N$  |

変数の右上につけて、どのエージェントとどのエージェントの取引もしくは貸借関係なのかを区別する

# 2. モデルの式とアルゴリズム
- $p^o=\frac{(1+\nu_1)(\sum_j W_{-1}^{j,o}+T_{v-1}^o+T_{c-1}^o+\delta k_{-1}^o)}{u^T\gamma_1 k_{-1}^o} + \nu_2(\lambda_5 \overline{p}_{-1}+\lambda_5 p_{-1}^o)$
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？
  - 収穫逓増の効果を入れるべき？入れるとしたらどうやって？
- 就業・失業判定、賃金率と賃金の決定のアルゴリズム
  - $EMP^j$ は失業時に0，就業時にoを値に持つ
  - if $\sum_o EMP_{-1}^{j}>0$
    - if $rand() < \zeta_1$
      - $EMP^{j}=0$
      - $w^j=w_{-1}^j[1 - \zeta_3 abs\{randn()\}]$
      - $W^{j,o} = v_{-1}^o w_{-1}^j$
    - else
      - $EMP^{j}=o=EMP_{-1}^{j}$
      - $w^j=w_{-1}^j[1 + \zeta_2 abs\{randn()\}]$
      - $W^{j,o} = v_{-1}^o w_{-1}^j$
  - else
    - 企業が Int64($\max[0, \frac{1}{A^o}\{u_{-1}^o k_{-1}^o-A_{-1}^o \sum_j (W_{-1}^{j,o}>0)\}]$) 人の求人を出す
    - 求人数に比例する確率で失業者は応募する。応募の中から定員までランダムに雇用する
    - $o'$で雇用が決まった場合
      - $w^j=w_{-1}^j[1 + \zeta_2 abs\{randn()\}]$
      - $EMP^{j}=o'$
    - 失業が続く場合
      - $w^j=w_{-1}^j[1 - \zeta_3 abs\{randn()\}]$
- $T_i^j=\tau_1 \sum_j (W_{-1}^{j,o}+P_{h-1}^{j,o}+S_{-1}^{j,n})$
- $T_a^j=\tau_2(\sum_o E_{h-1}^{o,j}+\sum_n M_{h-1}^{n,j}-\sum_n L_{h-1}^{j,n})$
- $T_v^o=\tau_3(\sum_j C_{-1}^{o,j}+I_{-1}^o+G_{-1}^o)$
- $T_c^o=\tau_4(\sum_j C_{-1}^{o,j}+G_{-1}^o+I_{-1}^o-\sum_j W_{-1}^{j,o}-T_v^o)$
- 政府支出のアルゴリズム
  - $G_{sum}=\sum G^o=G_0(1+\beta_1)^{(t-1)}$
  - $x^o=\max[1, x_{-1}^o \{1-\beta_2+\beta_3 randn\}]$
  - $z^o=\max(0, x^o-\beta_4)$
  - $G^o=\frac{z^o G_{sum}}{\sum z^o}$
  - 下限付きべき分布みたいなものを作って(xがこれにあたる)、一定の値を差し引いた後(zがこれに当たる)に、政府支出総額に合わせて大きくする。
- $g^o=\frac{G^o}{p^o}$
- $\sum_o C^{o,j}=\alpha_1(\sum_o W_{-1}^{o,j}-T_a^j-T_i^j-r_L \sum_n L_{h-1}^{j,n}+\sum_o P_{h-1}^{o,j}+\sum_n S_{-1}^{n,j}) + \alpha_2(\sum_o E_{h-1}^{o,j}+\sum_n M_{h-1}^{n,j}-\sum_n L_{h-1}^{n,j})$
- 消費財の購入先のアルゴリズム
  - 前回消費財を買った企業$o$から別の企業に乗り換える確率は
    - $rand() < \alpha_3 \frac{p^{o}-\overline{p}}{\overline{p}}$
    - とする。乗り換え先は「前期の消費財生産量($\sum_j c_{-1}^{o,j}$)＋前期の総消費量の一定割合($\alpha_4 \frac{\sum_o \sum_j c_{-1}^{o,j}+\sum_o g_{-1}^o}{J}$)」に比例する確率でランダムに決める
  - 今期の消費財購入先の企業を$o'$と書くことにすると、今期の消費水準は、
    - $c^{o'j}=\frac{\sum_o C^{o,j}}{p^{o'}}$
    - で決める。
    - $o'$以外のすべての企業の商品は消費しない($c^{o,j}=0 \ \ \ \{ o \mid o \neq o'\}$)
- $C^{o,j}=p^o c^{o,j}$
- $A^o=A_{-1}^o(1 + \mu_1 + \mu_2 \frac{i_{-1}^o}{k_{-1}^o})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u^o = \frac{\sum_j c^{o,j}+i^o+g^o}{\gamma_1 k_{-1}^o}$
- $i^o=\max\{0, \delta k_{-1}^o + (u^o-u^T)\gamma_1 k_{-1}^o + \gamma_2\frac{\sum_n M_{f-1}^{n,o}-\sum_n L_{f-1}^{o,n}}{p^o}\}$
- $I^o=p^o i^o$
- $k^o=(1-\delta)k_{-1}^o+i^o$
- $K^o=p^o k^o$
- $P^o=\sum_j C^{o,j}+G^o+I^o-\sum_j W^{j,o}-T_c^o-T_v^o-r_L \sum_n L_{f-1}^{o,n}$
- $P_h^{j,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{h-1}^{o,j}}{e_{-1}^o}$
- $P_b^{n,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{b-1}^{o,n}}{e_{-1}^o}$
- $P_f^o=P^o-\sum_j P_h^{j,o}-\sum_n P_b^{n,o}$
- $S^{j,n}=\{\theta_1(r_L L_{-1}^n+\sum_o P_b^{n,o})+\theta_2 \sum_o E_{b-1}^{o,n}\}\frac{f_{h-1}^{j,n}}{f_{-1}^n}$
- $NL_h^j=-\sum_o C^{o,j}+\sum_o W^{j,o}-T_i^j-T_a^j-r_L \sum_n L_{h-1}^{j,n}+\sum_o P_h^{j,o}+\sum_n S^{j,n}$
- $NL_f^o=-I^o+P_f^o$
- $NL_b^n=r_L L_{-1}^n+\sum_o P_b^{o,n}-\sum_j S^{j,n}$
- $NL_g=-\sum_o G^o+\sum_j T_i^j+\sum_j T_a^j+\sum_o T_v^o+\sum_o T_c^o$
- $\sum_o L_h^{j,n}=\epsilon_1 NL_h^j+\epsilon_2 \sum_o C^{o,j}$
  - 家計は預金口座を持つ銀行を借入先に選ぶ
    - $L_h^{j,n}=$
    - 本当は銀行ごとに信用スコアを計算して金利が安いところから借り入れようとするとか、既存の借入先を維持する傾向を持つとか、借入先を2つ以上にすることもあるとか、そういう効果を入れたいが、モデルの複雑さを抑えるためのアドホックな仮定として導入することにする。
- $\Delta L_h^{j,n}=L_h^{j,n}-L_{h-1}^{j,n}$
- $\sum_n \Delta L_f^{o,n}=\max\{-\sum_n L_{f-1}^{o,n}, (\lambda_3 + \lambda_4(\frac{P^o - P_f^o}{\sum_j E_{h-1}^{o,j} + \sum_n E_{b-1}^{o,n}} - r_L))(I^o+\sum_j W^{j,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{j,n} - \phi \sum_n M_{f-1}^{n,j})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
  - $L_h^{j,n}$ と同じく、企業が預金口座を持つ銀行から借り入れる
- $L_f^{o,n}=L_{f-1}^{o,n}+\Delta L_f^{o,n}$
- $L^n=\sum_j L_h^{j,n}+\sum_o L_f^{o,n}$
- $\Delta L^n=L^n-L_{-1}^n$
- $\Delta e^o=\frac{1}{p_{-1}^o}\{1-\lambda_3 - \lambda_4(\frac{P^o - P_f^o}{\sum_j E_{h-1}^{o,j} + \sum_n E_{b-1}^{o,n}} - r_L)\}(I^o+\sum_j W^{j,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{o,n} - \phi \sum_n M_{f-1}^{n,o})$
- $E^o=E_{-1}^o+p_{e-1}^o\Delta e^o$
- $e^o=e_{-1}^o+\Delta e^o$
- 家計はどうやって株式保有額を決めるのか
  - どこの株を保有するかの意思決定
    - 利益、純金融資産、配当率　で重み付けされたランダムサンプリングと、前期保有していた銘柄を持ち続ける傾向の合成
  - $E_h^j=$
  - $F_h^j=$
- 銀行の企業株式保有割合を決める
  - 利益、純金融資産、配当率　で重み付けされた保有割合に、摂動を加え、期首の保有割合の影響も残す
  - $E_b^{o,n}=$
- $p_e^o=\frac{\sum_j E_h^{o,j}+\sum_n E_b^{o,n}}{e^o}$
- $e_h^{o,j}=\frac{E_h^{o,j}}{p_e^o}$
- $e_b^{o,n}=\frac{E_b^{o,n}}{p_e^o}$
- $\Delta e_h^{o,j}=e_h^{o,j}-e_{h-1}^{o,j}$
- $\Delta e_b^{o,n}=e_h^{o,n}-e_{b-1}^{o,n}$
- 家計の預金の決定
  - $\sum_n \Delta M_h^{n,j}=NL_h^j-\sum_o(p_{e-1}^o\Delta e_h^{o,j})+\sum_n \Delta L_h^{j,n}$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_h=M_{h-1}+\Delta M_h$
- 企業の預金の決定
  - $\sum_n \Delta M_f^{n,o}=NL_f^o+\sum_n \Delta L_f^{o,n}+p_{e-1}^o\Delta e^o$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_f^{n,o}=M_{f-1}^{n,o}+\Delta M_f^{n,o}$
- $M^n=\sum_j M_h^{n,j}+\sum_o M_f^{n,o}$
- $\Delta M^n=M^n-M_{-1}^n$
- $\Delta H^n = NL_b^n - \sum_o p_{e-1}^o\Delta e_b^{o,n} + \Delta M^n - \Delta L^n$
- $H^n = H_{-1}^n + \Delta H^n$
- $NW_h^j=\sum_n M_h^{n,j}-\sum_n L_h^{j,n}+\sum_o E_h^{o,j}$
- $NW_f^o=K^o+\sum_n M_f^{n,o}-\sum_n L_f^{o,n}-E^o$
- $NW_b^n=-M^n+L^n+\sum_o E_b^{o,n}+H^n$
- $NW_g=-\sum_n H^n$
- $D_E=-\sum_o E^o+\sum_o\sum_j E_h^{o,j}+\sum_o\sum_j E_b^{o,j}$
- $D_F=\sum_j\sum_n F_h^{j,n}-\sum_n F^n$
- $NW=\sum_j NW_h^j+\sum_o NW_f+\sum_n NW_b+NW_g$
- $v^o=\frac{u^o k^o}{A^o \sum_j (W^{j,o}>0)}$
- 
- 倒産と新規参入のアルゴリズムも追加する必要がある
- 倒産
  - 純金融資産が0以下になると倒産する
  - 

隠された恒等式として

- 

を使ってモデルの会計的一貫性を確認する


# 拡張の方向性
## アルゴリズム
- 3次元配列を用いず、借方に現れるエージェント番号xをインデックスとする二次元配列を定義(arr[x,2])し、二次元目の1つ目の要素に借方に現れるエージェント番号を、2つ目の要素に金額を、記録する方法を検討する
  - メモリ消費量の節約のため
## 行動の仮定
- 政府や銀行が家計を雇う
- 国債の導入
- 銀行への資本規制の追加
- 金利の内生化
- 預金/借入金の多対多化
- 貸し渋り、貸し剥がしのような与信行動を追加する
- 失業者への社会保障をモデル含める
- など