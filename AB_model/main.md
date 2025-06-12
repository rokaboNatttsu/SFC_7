- [1. エージェント番号](#1-エージェント番号)
- [2. モデルの式とアルゴリズム](#2-モデルの式とアルゴリズム)
- [拡張の方向性](#拡張の方向性)



AB_simulatio,n.jlで実装する予定

# 1. エージェント番号
|      | インデックス | 総数 |
| :--- | :----------: | :--: |
| 家計 |     $l$      | $L$  |
| 企業 |     $o$      | $O$  |
| 銀行 |     $n$      | $N$  |

変数の右上につけて、どのエージェントとどのエージェントの取引もしくは貸借関係なのかを区別する

# 2. モデルの式とアルゴリズム
- $p^o=\lambda_p \frac{(1+\nu_1+\nu_2\frac{\sum_n L_{f-1}^{o,n}}{\sum_l C_{-1}^{o,l}+ G_{-1}^o})(\sum_l W_{-1}^{l,o}+T_{v-1}^o+T_{c-1}^o+\delta k_{-1}^o)}{u^T\gamma_1 k_{-1}^o} + (1-\lambda_p)\nu_3(p_{-1}^o-\overline{p}_{-1})$
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？
  - 価格競争の成分は？ $\overline{p}_{-1}=\frac{\sum_o(\sum_l C^{o,l}_{-1}+G_{-1}^o)}{\sum_o(\sum_l c^{o,l}_{-1}+g_{-1}^o)}$ に近づくように $\nu_3(p_{-1}^o-\overline{p}_{-1})$ の項を追加する？
  - 収穫逓増の効果を入れるべき？入れるとしたらどうやって？
  - $\lambda_p$を適応ラメータにして、消費者が価格差に対しする商品の乗り換え速度に適応するようにしたほうが良いかも
- $T_i^l=\tau_1 \sum_l (W_{-1}^{l,o}+P_{h-1}^{l,o}+S_{-1}^{l,n})$
- $T_a^l=\tau_2(\sum_o E_{h-1}^{o,l}+\sum_n M_{h-1}^{n,l}-\sum_n L_{h-1}^{l,n})$
- $T_v^o=\tau_3(\sum_l C_{-1}^{o,l}+I_{-1}^o+G_{-1}^o)$
- $T_c^o=\tau_4(\sum_l C_{-1}^{o,l}+G_{-1}^o+I_{-1}^o-\sum_l W_{-1}^{l,o}-T_v^o)$
- 政府支出のアルゴリズム　govExpAlg2.jl のアルゴリズムを採用
  - $G_{sum}=\sum G^o=G_0(1+\beta_1)^{(t-1)}$
  - $x^o=\max[1, x_{-1}^o \{1+\beta_2(-1+\beta_3 randn)\}]$
  - $z^o=\max(0, x^o-\beta_4)$
  - $G^o=\frac{z^o G_{sum}}{\sum z^o}$
  - 下限付きべき分布みたいなものを作って(xがこれにあたる)、一定の値を差し引いた後(zがこれに当たる)に、政府支出総額に合わせて大きくする。
- $g^o=\frac{G^o}{p^o}$
- $\sum_o C^{o,l}=\alpha_1(\sum_o W_{-1}^{o,l}-T_a^l-T_i^l-r_L \sum_n L_{h-1}^{l,n}+\sum_o P_{h-1}^{o,l}+\sum_n S_{-1}^{n,l}) + \alpha_2(\sum_o E_{h-1}^{o,l}+\sum_n M_{h-1}^{n,l}-\sum_n L_{h-1}^{n,l})$
- 消費財の購入先のアルゴリズム
  - 前回消費財を買った企業$o$から別の企業に乗り換える確率は
    - $rand() < \alpha_3 \frac{p^{o}-\overline{p}}{\overline{p}}$
    - とする。乗り換え先は「前期の消費財生産量($\sum_l c_{-1}^{o,l}$)＋前期の総消費量の一定割合($\alpha_4 \sum_o \sum_l c_{-1}^{o,l}$)」に比例する確率でランダムに決める
  - 今期の消費財購入先の企業を$o'$と書くことにすると、今期の消費水準は、
    - $c^{o'l}=\frac{\sum_o C^{o,l}}{p^{o'}}$
    - で決める。
    - $o'$以外のすべての企業の商品は消費しない($c^{o,l}=0 \ \ \ \{ o \mid o \neq o'\}$)
- $C^{o,l}=p^o c^{o,l}$
- $A^o=A_{-1}^o(1 + \mu_1 + \mu_2 \frac{i_{-1}^o}{k_{-1}^o})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u^o = \frac{\sum_l c^{o,l}+i^o+g^o}{\gamma_1 k_{-1}^o}$
- $i^o=\delta k_{-1}^o + (u_{-1}^o-u^T)\gamma_2 k_{-1}^o + \gamma_3\frac{\sum_n M_{f-1}^{n,o}-\sum_n L_{f-1}^{o,n}}{p^o}$
- $I^o=p^o i^o$
- $k^o=(1-\delta)k_{-1}^o+i^o$
- $K^o=p^o k^o$
- 就業・失業判定アルゴリズム　外部にファイルを作る
- $EMP^l$ は失業時に0，就業時にoを値に持つ
  - if $\sum_o EMP_{-1}^{l}>0$
    - if $rand() < \zeta_3$
      - $EMP^{l}=0$
      - 失業者への社会保障をモデル含めることを検討
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}^l[1 - \zeta_3 abs\{randn()\}]$
    - else
      - $EMP^{l}=EMP_{-1}^{l}$
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}[1 + \zeta_2 abs\{randn()\}]$
      - $W^{l,EMP^l} = w^l$
  - else
    - 企業が Int64($\max[0, \frac{1}{A^o}\{u_{-1}^o k_{-1}^o-A^o \sum_l (w_{-1}^{l,o}>0)\}]$) 人の求人を出す
    - 求人数に比例する確率で失業者は応募する。応募の中から定員までランダムに雇用する
    - $o'$で雇用が決まった場合
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}^l[1 + \zeta_2 abs\{randn()\}]$
      - $EMP^{l}=o'$
      - $W^{l,EMP^l} = w^l$
    - 失業が続く場合
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}^l[1 - \zeta_3 abs\{randn()\}]$
- $P^o=\sum_l C^{o,l}+G^o+I^o-\sum_l W^{l,o}-T_c^o-T_v^o-r_L \sum_n L_{f-1}^{o,n}$
- $P_h^{l,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{h-1}^{o,l}}{e_{-1}^o}$
- $P_b^{n,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{b-1}^{o,n}}{e_{-1}^o}$
- $P_f^o=P^o-\sum_l P_h^{l,o}-\sum_n P_b^{n,o}$
- $S^{l,n}=\{\theta_3(r_L L_{-1}^n+\sum_o P_b^{n,o})+\theta_4 \sum_o E_{b-1}^{o,n}\}\frac{f_{h-1}^{l,n}}{f_{-1}^n}$
- $NL_h^l=-\sum_o C^{o,l}+\sum_o W^{l,o}-T_i^l-T_a^l-r_L \sum_n L_{h-1}^{l,n}+\sum_o P_h^{l,o}+\sum_n S^{l,n}$
- $NL_f^o=-I^o+P_f^o$
- $NL_b^n=r_L L_{-1}^n+\sum_o P_b^{o,n}-\sum_l S^{l,n}$
- $NL_g=-\sum_o G^o+\sum_l T_i^l+\sum_l T_a^l+\sum_o T_v^o+\sum_o T_c^o$
- $\sum_o L_h^{l,n}=\epsilon_1 NL_h^l+\epsilon_2 \sum_o C^{o,l}$
  - 家計は預金口座を持つ銀行を借入先に選ぶ
    - $L_h^{l,n}=$
    - 本当は銀行ごとに信用スコアを計算して金利が安いところから借り入れようとするとか、既存の借入先を維持する傾向を持つとか、借入先を2つ以上にすることもあるとか、そういう効果を入れたいが、モデルの複雑さを抑えるためのアドホックな仮定として導入することにする。
- $\Delta L_h^{l,n}=L_h^{l,n}-L_{h-1}^{l,n}$
- $\sum_n \Delta L_f^{o,n}=\max\{-\sum_n L_{f-1}^{o,n}, (\lambda_3 + \lambda_4(\frac{P^o - P_f^o}{\sum_l E_{h-1}^{o,l} + \sum_n E_{b-1}^{o,n}} - r_L))(I^o+\sum_l W^{l,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{l,n} - \phi \sum_n M_{f-1}^{n,l})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
  - $L_h^{l,n}$ と同じく、企業が預金口座を持つ銀行から借り入れる
- $L_f^{o,n}=L_{f-1}^{o,n}+\Delta L_f^{o,n}$
- $L^n=\sum_l L_h^{l,n}+\sum_o L_f^{o,n}$
- $\Delta L^n=L^n-L_{-1}^n$
- $\Delta e^o=\frac{1}{p_{-1}^o}\{1-\lambda_3 - \lambda_4(\frac{P^o - P_f^o}{\sum_l E_{h-1}^{o,l} + \sum_n E_{b-1}^{o,n}} - r_L)\}(I^o+\sum_l W^{l,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{o,n} - \phi \sum_n M_{f-1}^{n,o})$
- $E^o=E_{-1}^o+p_{e-1}^o\Delta e^o$
- $e^o=e_{-1}^o+\Delta e^o$
- 家計はどうやって株式保有額を決めるのか
  - どこの株を保有するかの意思決定
    - 利益、純金融資産、配当率　で重み付けされたランダムサンプリングと、前期保有していた銘柄を持ち続ける傾向の合成
  - $E_h^l=$
  - $F_h^l=$
- 銀行の企業株式保有割合を決める
  - 利益、純金融資産、配当率　で重み付けされた保有割合に、摂動を加え、期首の保有割合の影響も残す
  - $E_b^{o,n}=$
- $p_e^o=\frac{\sum_l E_h^{o,l}+\sum_n E_b^{o,n}}{e^o}$
- $e_h^{o,l}=\frac{E_h^{o,l}}{p_e^o}$
- $e_b^{o,n}=\frac{E_b^{o,n}}{p_e^o}$
- $\Delta e_h^{o,l}=e_h^{o,l}-e_{h-1}^{o,l}$
- $\Delta e_b^{o,n}=e_h^{o,n}-e_{b-1}^{o,n}$
- 家計の預金の決定
  - $\sum_n \Delta M_h^{n,l}=NL_h^l-\sum_o(p_{e-1}^o\Delta e_h^{o,l})+\sum_n \Delta L_h^{l,n}$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_h=M_{h-1}+\Delta M_h$
- 企業の預金の決定
  - $\sum_n \Delta M_f^{n,o}=NL_f^o+\sum_n \Delta L_f^{o,n}+p_{e-1}^o\Delta e^o$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_f^{n,o}=M_{f-1}^{n,o}+\Delta M_f^{n,o}$
- $M^n=\sum_l M_h^{n,l}+\sum_o M_f^{n,o}$
- $\Delta M^n=M^n-M_{-1}^n$
- 
- 
- ここからAB化の作業を再開する
- $\Delta H = NL_b - p_{e-1}\Delta e_b + \Delta M - \Delta L$
- $H = H_{-1} + \Delta H$
- $NW_h=M_h-L_h+E_h$
- $NW_f=K+M_f-L_f-E$
- $NW_b=-M+L+E_b+H$
- $NW_g=H$
- $D=-E+E_h+E_b$
- $NW=NW_h+NW_f+NW_b+NW_g$
- 倒産と新規参入のアルゴリズムも追加する必要がある

隠された恒等式として

- $K+D-NW=0$
- $NL_h+NL_f+NL_b+NL_g=0$
- $\Delta L=\Delta L_h+\Delta L_f$
- $\Delta M=\Delta M_h+\Delta M_f$
- $NL_g+\Delta H=0$

を使ってモデルの会計的一貫性を確認する


# 拡張の方向性
- 政府や銀行が家計を雇う
- 国債の導入
- 銀行への資本規制の追加
- 金利の内生化
- 預金/借入金の多対多化
- など