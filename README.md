- [1. 集計的モデル](#1-集計的モデル)
  - [1.1. TFM](#11-tfm)
  - [1.2. BSM](#12-bsm)
  - [1.3. モデルの式一覧](#13-モデルの式一覧)
- [2. エージェント・ベース・モデル](#2-エージェントベースモデル)
  - [2.1. モデルの説明](#21-モデルの説明)

# 1. 集計的モデル

いずれ作るAB-SFCモデルの基礎とする予定

simulation.ipynbで実装した

拡張する方向として検討すること

- 公債と公債利回りの追加
  - 公債利回りと民間への貸出金利の適応的連動 $r_L=(1-\lambda_r)r_{L-1}+\lambda_r \epsilon r_{GB}$ 
- 労働者と投資家の分割
- 公務員の追加
- 政府の資本の追加
- 資本財生産企業と消費財生産企業の分割
- 海外の導入

## 1.1. TFM
|              |         家計          |   企業(経常)   |     企業(資本)      |         銀行          |  統合政府   |  合計   |
| :----------- | :-------------------: | :------------: | :-----------------: | :-------------------: | :---------: | :-----: |
| 消費         |         $-C$          |      $+C$      |                     |                       |             |   $0$   |
| 政府支出     |                       |      $+G$      |                     |                       |    $-G$     |   $0$   |
| 投資         |                       |      $+I$      |        $-I$         |                       |             |   $0$   |
| 賃金         |         $+W$          |      $-W$      |                     |                       |             |   $0$   |
| 所得税       |        $-T_i$         |                |                     |                       |   $+T_i$    |   $0$   |
| 資産税       |        $-T_a$         |                |                     |                       |   $+T_a$    |   $0$   |
| 付加価値税   |                       |     $-T_v$     |                     |                       |   $+T_v$    |   $0$   |
| 法人税       |                       |     $-T_c$     |                     |                       |   $+T_c$    |   $0$   |
| 借入金利払   |    $-r_L L_{h-1}$     | $-r_L L_{f-1}$ |                     |     $+r_L L_{-1}$     |             |   $0$   |
| 企業利潤     |        $+P_h$         |      $-P$      |       $+P_f$        |        $+P_b$         |             |   $0$   |
| 銀行配当     |         $+S$          |                |                     |         $-S$          |             |   $0$   |
| [ 純貸出 ]   |      [ $NL_h$ ]       |                |     [ $NL_f$ ]      |      [ $NL_b$ ]       | [ $NL_g$ ]  | [ $0$ ] |
| 株式純発行   | $-p_{e-1} \Delta e_h$ |                | $+p_{e-1} \Delta e$ | $-p_{e-1} \Delta e_b$ |             |   $0$   |
| 預金の移動   |     $-\Delta M_h$     |                |    $-\Delta M_f$    |      $+\Delta M$      |             |   $0$   |
| 貸付金の移動 |     $+\Delta L_h$     |                |    $+\Delta L_f$    |      $-\Delta L$      |             |   $0$   |
| 現金の移動   |                       |                |                     |      $-\Delta H$      | $+\Delta H$ |   $0$   |
| 合計         |          $0$          |      $0$       |         $0$         |          $0$          |     $0$     |         |

## 1.2. BSM
|        |  家計   |  企業   |  銀行   | 統合政府 | 合計  |
| :----- | :-----: | :-----: | :-----: | :------: | :---: |
| 資本   |         |  $+K$   |         |          | $+K$  |
| 預金   | $+M_h$  | $+M_f$  |  $-M$   |          |  $0$  |
| 貸付金 | $-L_h$  | $-L_f$  |  $+L$   |          |  $0$  |
| 株式   | $+E_h$  |  $-E$   | $+E_b$  |          | $+D$  |
| 現金   |         |         |  $+H$   |   $-H$   |  $0$  |
| 純資産 | $-NW_h$ | $-NW_f$ | $-NW_b$ | $-NW_g$  | $-NW$ |
| 合計   |   $0$   |   $0$   |   $0$   |   $0$    |  $0$  |

## 1.3. モデルの式一覧
- $p=\frac{(1+\nu_1+\nu_2\frac{L_{f-1}}{C_{-1}+G_{-1}})(W_{-1}+T_{v-1}+T_{c-1})}{u^T\gamma_1 k_{-1}}$
  - 価格を安定させたがること（塩沢）を反映した価格設定式はないか？
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？集計的SFCモデルでは投資資金回収率に比例する項と定数項の和が価格マークアップの水準になるのが現実的か？
- $T_i=\tau_1 W_{-1}$
- $T_a=\tau_2(E_{h-1}+M_{h-1}-L_{h-1})$
- $T_v=\tau_3(C_{-1}+I_{-1}+G_{-1})$
- $T_c=\tau_4(C_{-1}+G_{-1}+I_{-1}-W_{-1}-T_v)$
- $G=(1+\beta)G_{-1}$
- $g=\frac{G}{p}$
- $C=\alpha_1(W_{-1}-T_a-T_i-r_L L_{-1}+P_{h-1}+S_{-1}) + \alpha_2(E_{h-1}+M_{-1}-L_{-1})$
- $c=\frac{C}{p}$
- $A=A_{-1}(1 + \mu_1 + \mu_2 \frac{i_{-1}}{k_{-1}})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u = \frac{c+i+g}{\gamma_1 k_{-1}}$
- $i=\delta k_{-1} + (u_{-1}-u^T)\gamma_2 k_{-1} + \gamma_3\frac{M_{f-1}-L_{f-1}}{p}$
- $I=p i$
- $k=(1-\delta)k_{-1}+i$
- $K=p k$
- $v=u \frac{k_{-1}}{A}$
  - 労働稼働率は、技術水準（と労働者数の積）/消費および投資の実質需要　で定義する。
- $w=(1-\zeta_1)w_{-1} + \zeta_1 w_{-1}\{1 + \zeta_2(v - v^T)\}$
- $W=w(c_{-1}+g_{-1}+i_{-1})$
  - 労働者側の賃金交渉力が高いときは名目賃金が上昇しやすく、逆に労働者側の賃金交渉力が低いときは名目賃金が上がりにくい。労働者側の賃金交渉力の相対的な強さは、労働稼働率の一次関数で決まることにする。
- $P=C+G+I-W-T_c-T_v-r_L L_{f-1}$
- $P_h=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{h-1}}{e_{-1}}$
- $P_b=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{b-1}}{e_{-1}}$
- $P_f=P-P_h-P_b$
- $S=\theta_3(r_L L_{-1}+P_b)+\theta_4 E_{b-1}$
- $NL_h=-C+W-T_i-T_a-r_L L_{h-1}+P_h+S$
- $NL_f=-I+P_f$
- $NL_b=r_L L_{-1}+P_b-S$
- $NL_g=-G+T_i+T_a+T_v+T_c$
- $L_h=\epsilon_1 NL_h+\epsilon_2 C$
- $\Delta L_h=L_h-L_{h-1}$
- $\Delta L_f=\max\{-L_{f-1}, (\lambda_3 + \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
- $L_f=L_{f-1}+\Delta Lf$
- $L=L_h+L_f$
- $\Delta L=L-L_{-1}$
- $\Delta e=\frac{1}{p_{-1}}(1-\lambda_3 - \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})$
- $E=E_{-1}+p_{e-1}\Delta e$
- $e=e_{-1}+\Delta e$
- $E_h=(\lambda_1+\lambda_2\frac{P_{h-1}}{E_{h-1}})(NL_h+\Delta L_h+E_{h-1}+M_{h-1})$
- $E_b=(1-\lambda_5 r_L+\lambda_5 \frac{P_{b-1}}{E_{b-1}})(NL_b+E_{b-1})$
- $p_e=\frac{E_h+E_b}{e}$
- $e_h=\frac{E_h}{p_e}$
- $e_b=\frac{E_b}{p_e}$
- $\Delta e_h=e_h-e_{h-1}$
- $\Delta e_b=e_h-e_{b-1}$
- $\Delta M_h=NL_h-p_{e-1}\Delta e_h+\Delta L_h$
- $M_h=M_{h-1}+\Delta M_h$
- $\Delta M_f=NL_f+\Delta L_f+p_{e-1}\Delta e$
- $M_f=M_{f-1}+\Delta M_f$
- $M=M_h+M_f$
- $\Delta M=M-M_{-1}$
- $\Delta H = NL_b - p_{e-1}\Delta e_b + \Delta M - \Delta L$
- $H = H_{-1} + \Delta H$
- $NW_h=M_h-L_h+E_h$
- $NW_f=K+M_f-L_f-E$
- $NW_b=-M+L+E_b+H$
- $NW_g=H$
- $D=-E+E_h+E_b$
- $NW=NW_h+NW_f+NW_b+NW_g$

隠された恒等式として

- $K+D-NW=0$
- $NL_h+NL_f+NL_b+NL_g=0$
- $\Delta L=\Delta L_h+\Delta L_f$
- $\Delta M=\Delta M_h+\Delta M_f$
- $NL_g+\Delta H=0$

を使ってモデルの会計的一貫性を確認する

# 2. エージェント・ベース・モデル
|      | インデックス | 総数 |
| :--- | :----------: | :--: |
| 家計 |     $l$      | $L$  |
| 企業 |     $m$      | $M$  |
| 銀行 |     $n$      | $N$  |

変数の右上につけて、エージェント番号を区別する

## 2.1. モデルの説明
- $p^w=\frac{(1+\nu_1+\nu_2\frac{L_{f-1}^m}{\sum_l C_{-1}^{ml}+ G_{-1}^m})(\sum_l W_{-1}^{ml}+T_{v-1}^m+T_{c-1}^m)}{u^T\gamma_1 k_{-1}^m}$
  - 価格を安定させたがること（塩沢）を反映した価格設定式はないか？
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？集計的SFCモデルでは投資資金回収率に比例する項と定数項の和が価格マークアップの水準になるのが現実的か？
- $T_i^l=\tau_1 \sum_l W_{-1}^{ml}$
- $T_a^l=\tau_2(\sum E_{h-1}^m+\sum M_{h-1}^n-\sum L_{h-1}^n)$
- $T_v^m=\tau_3(\sum_l C_{-1}^{lm}+I_{-1}^m+G_{-1}^m)$
- $T_c=\tau_4(\sum_l C_{-1}^{lm}+G_{-1}^m+I_{-1}^m-\sum_l W_{-1}^{lm}-T_v^m)$
- 政府支出のアルゴリズム
  - 公共事業発注先を指定した数入れ替える。新たな契約を結ぶときは、その金額はその期の政府支出の一定割合にする（増加率が歩調を揃える）
  - 発注金額を各企業ごとにランダムに上下させる。
  - $\sum G^m=(1+\beta)\sum G_{-1}^m$ の束縛条件を満たすよう、発注金を一定倍にして調整する
  - govExpAlg.jl でアルゴリズムを書いた
  - 規模が一極集中し安すぎる問題。特定のべき分布（下限あり）に収束するような規則を考え出すか、[対数正規分布に従う移動規則](https://rokabonatttsu.hatenablog.com/entry/2023/01/25/215420)を使うことを検討する
  - 大きな発注先との取引をやめると、他の企業の稼働率が上がり、大量失業がそれに続くことになるのでは？という懸念がある
- 
- $g=\frac{G}{p}$
- $C=\alpha_1(W_{-1}-T_a-T_i-r_L L_{-1}+P_{h-1}+S_{-1}) + \alpha_2(E_{h-1}+M_{-1}-L_{-1})$
- $c=\frac{C}{p}$
- $A=A_{-1}(1 + \mu_1 + \mu_2 \frac{i_{-1}}{k_{-1}})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u = \frac{c+i+g}{\gamma_1 k_{-1}}$
- $i=\delta k_{-1} + (u_{-1}-u^T)\gamma_2 k_{-1} + \gamma_3\frac{M_{f-1}-L_{f-1}}{p}$
- $I=p i$
- $k=(1-\delta)k_{-1}+i$
- $K=p k$
- $v=u \frac{k_{-1}}{A}$
  - 労働稼働率は、技術水準（と労働者数の積）/消費および投資の実質需要　で定義する。
- $w=(1-\zeta_1)w_{-1} + \zeta_1 w_{-1}\{1 + \zeta_2(v - v^T)\}$
- $W=w(c_{-1}+g_{-1}+i_{-1})$
  - 労働者側の賃金交渉力が高いときは名目賃金が上昇しやすく、逆に労働者側の賃金交渉力が低いときは名目賃金が上がりにくい。労働者側の賃金交渉力の相対的な強さは、労働稼働率の一次関数で決まることにする。
- $P=C+G+I-W-T_c-T_v-r_L L_{f-1}$
- $P_h=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{h-1}}{e_{-1}}$
- $P_b=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{b-1}}{e_{-1}}$
- $P_f=P-P_h-P_b$
- $S=\theta_3(r_L L_{-1}+P_b)+\theta_4 E_{b-1}$
- $NL_h=-C+W-T_i-T_a-r_L L_{h-1}+P_h+S$
- $NL_f=-I+P_f$
- $NL_b=r_L L_{-1}+P_b-S$
- $NL_g=-G+T_i+T_a+T_v+T_c$
- $L_h=\epsilon_1 NL_h+\epsilon_2 C$
- $\Delta L_h=L_h-L_{h-1}$
- $\Delta L_f=\max\{-L_{f-1}, (\lambda_3 + \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
- $L_f=L_{f-1}+\Delta Lf$
- $L=L_h+L_f$
- $\Delta L=L-L_{-1}$
- $\Delta e=\frac{1}{p_{-1}}(1-\lambda_3 - \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})$
- $E=E_{-1}+p_{e-1}\Delta e$
- $e=e_{-1}+\Delta e$
- $E_h=(\lambda_1+\lambda_2\frac{P_{h-1}}{E_{h-1}})(NL_h+\Delta L_h+E_{h-1}+M_{h-1})$
- $E_b=(1-\lambda_5 r_L+\lambda_5 \frac{P_{b-1}}{E_{b-1}})(NL_b+E_{b-1})$
- $p_e=\frac{E_h+E_b}{e}$
- $e_h=\frac{E_h}{p_e}$
- $e_b=\frac{E_b}{p_e}$
- $\Delta e_h=e_h-e_{h-1}$
- $\Delta e_b=e_h-e_{b-1}$
- $\Delta M_h=NL_h-p_{e-1}\Delta e_h+\Delta L_h$
- $M_h=M_{h-1}+\Delta M_h$
- $\Delta M_f=NL_f+\Delta L_f+p_{e-1}\Delta e$
- $M_f=M_{f-1}+\Delta M_f$
- $M=M_h+M_f$
- $\Delta M=M-M_{-1}$
- $\Delta H = NL_b - p_{e-1}\Delta e_b + \Delta M - \Delta L$
- $H = H_{-1} + \Delta H$
- $NW_h=M_h-L_h+E_h$
- $NW_f=K+M_f-L_f-E$
- $NW_b=-M+L+E_b+H$
- $NW_g=H$
- $D=-E+E_h+E_b$
- $NW=NW_h+NW_f+NW_b+NW_g$

隠された恒等式として

- $K+D-NW=0$
- $NL_h+NL_f+NL_b+NL_g=0$
- $\Delta L=\Delta L_h+\Delta L_f$
- $\Delta M=\Delta M_h+\Delta M_f$
- $NL_g+\Delta H=0$

を使ってモデルの会計的一貫性を確認する
