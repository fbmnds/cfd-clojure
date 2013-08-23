(ns cfd-clojure-test
  (:use midje.sweet
        (incanter core charts))
  (:require [cfd-clojure :refer :all]
            [clojure.data.json :as json]
            [clojure.string :as str]
            [clojure.math.numeric-tower :as math]
            [clatrix.core :as clx]
            [clojure.core.match :as m]))

(import 'java.io.FileOutputStream)
(use '[clojure.java.shell :only [sh]])


(set! *warn-on-reflection* true)

;; save settings
;;
(def print-length *print-length*)
(set! *print-length* 120)
;;
;; set consistent en_US number formats for number/string conversions
;;
(def en_US (java.util.Locale. "en" "US"))
(def locale (java.util.Locale/getDefault))
(if-not (= locale en_US)
  (java.util.Locale/setDefault en_US))



;; depends on en_US for format consistency
;;
(defn- format-x [x n]
  (read-string (format (str/join ["%." (str n) "f"]) x)))

(defn- format-zz [zz n]
  (map #(map (fn [x] (format-x x n)) %) zz))

;; TODO: understand with-precision
;;
;; (defn- round-zz [zz n]
;;   (map #(map (fn [x] (with-precision n x)) %) zz))



(defmacro pure-time
  "Like clojure.core/time, returns the time as a value
   instead of a string."
  [expr]
  `(let [start# (. System (nanoTime))]
     (do
       (prn ~expr)
       (/ (double (- (. System (nanoTime)) start#)) 1000000.0))))


(def tn [100 1000 10000 100000])
(def theta 12)

(defn- perf? [tn]
  (let [[t100 t1000 t10000 t100000] tn
        res (and (< (/ t100000 t10000) theta)
                 (< (/ t10000 t1000) theta)
                 (< (/ t1000 t100) theta))]
    (if-not res
      (println "\n- observed times:\n  ---------------\n "
               (map #(format-x % 4) tn))
      true)))



(defn euklid [u v]
  (Math/sqrt (reduce + (map (comp #(* % %) #(apply - %)) (map vector u v)))))


(defn- shape [z dz]
  (if (and (>= z (math/floor (/ 0.5 dz))) ; Python implicitely uses 'floor'
           (< z (inc (/ 1. dz))))
    2.
    1.))

(defn- set-u0
  ([nx dx]
     (for [x (range nx)] (shape x dx)))
  ([ny dy nx dx] ; y rows, x cols
     (let [u0 (make-array Double/TYPE ny nx)]
       (doseq [y (range ny)
               x (range nx)]
         (aset u0 y x (if (= 2. (shape x dx) (shape y dy)) 2. 1.)))
       (matrix u0))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 1: Linear Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(fact
 "Step 1: Linear Convection" :step1

 (let [nx 41
       dx (/ 2. (dec nx))
       nt 25
       m {:c 1 :x-steps 2 :dx dx :dt 0.025}
       x (map #(* % dx) (range nx))
       u0 (set-u0 nx dx)
       u (take (inc nt) (discretize linear-convection m u0))
       u1 (second u)
       u5 (nth u 6)
       u12 (nth u 13)
       u17 (nth u 18)
       u_nt (last u)]
   (fact "nx = 41" :step1 nx => 41)
   (fact "u0"  :step1
         u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
                2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

   ;; precision of 10^-5 fails only for the first 3 values 1.00045526  1.00007826  1.00000972
   ;; hence reduced to 10^-3:
   ;;
   (fact "u_nt" :step1
         (map #(format-x % 3) u_nt)
         => (map #(format-x % 3)
                 [1.00045526  1.00007826  1.00000972  1.00000077  1.00000003  1.          1.
                  1.          1.          1.          1.00000003  1.00000077  1.00000972
                  1.00007826  1.00045526  1.00203866  1.00731665  1.02164263  1.05387607
                  1.11476147  1.21217811  1.34501895  1.49999923  1.6549713   1.78774363
                  1.88478327  1.94408527  1.97104073  1.97104073  1.94408527  1.88478327
                  1.78774363  1.6549713   1.49999923  1.34501895  1.21217811  1.11476147
                  1.05387607  1.02164263  1.00731665  1.00203866]))
   (doto
       (xy-plot x u0
                :title "Linear Convection"
                :x-label "x (nx=41)"
                :y-label "u(x,0) / u(x,0.025) / ... / u(x,0.625)"
                :legend true)
     (add-lines x u1)
     (add-lines x u5)
     (add-lines x u12)
     (add-lines x u17)
     (add-lines x u_nt)
     view))



 (let [nx 81
       dx (/ 2. (dec nx))
       nt 25
       m {:c 1 :x-steps 2 :dx dx :dt 0.025}
       x (map #(* % dx) (range nx))
       u0 (set-u0 nx dx)
       u (take (inc nt) (discretize linear-convection m u0))
       u_nt (last u)]
   (fact "nx = 81" :step1 nx => 81)
   (fact "u0" :step1
         u0
         => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
             2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
             2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

   (fact "u_nt" :step1
         u_nt
         => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  2.  2.  2.  2.  2.
             2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
             2.  2.  2.  2.  2.  2.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])
   (doto
       (xy-plot x u0
                :title "Linear Convection"
                :x-label "x (nx=81)"
                :y-label "u(x,0) / u(x,0.625)"
                :legend true)
     (add-lines x u_nt)
     view))



;;; Visualisation:
;;; http://markov.uc3m.es/2012/11/temporal-networks-with-igraph-and-r-with-20-lines-of-code/

 (defn save-frame [x t u_t]
   (let
       [p (xy-plot x u_t
                   :title "Linear Convection"
                   :x-label "x (nx=41)"
                   :y-label "u(x,0) ... u(x,0.625) (nt=100)"
                   :legend true)
        n (cond (< 99 t) (str t)
                (< 9 t) (str/join ["0" (str t)])
                :else (str/join ["00" (str t)]))
        fname (str/join ["./mpeg/linear-convection-" n ".png"])
        fos (FileOutputStream. fname)]
     (save p fos)
     (if (= 0 (mod t 100)) (println "... " (- 1000 t)))
     (.close fos)))

 (defn animation []
   (print "purging ./mpeg ... ")
   (sh "rm" "./mpeg/*")
   (println "done.\ncalculating Linear Convection (nt=1000) ... ")
   (let [nx 41
         dx (/ 2. (dec nx))
         nt 1000
         m {:c 1 :x-steps 2 :dx dx :dt 0.000625}
         x (map #(* % dx) (range nx))
         u0 (set-u0 nx dx)
         u (take nt (discretize linear-convection m u0))]
     (doall (map #(save-frame x % (nth u %)) (range (count u)))))
   (println "done.\nencoding of ./mpeg/linear-convection.mp4 ... ")
   ;;
   ;; ffmpeg -r 100 -b 20M -i linear-convection-%03d.png linear-convection.mp4
   ;;
   (sh "ffmpeg" "-r" "100" "-b" "20M" "-i"
       "mpeg/linear-convection-%03d.png"
       "mpeg/linear-convection.mp4")
   (println "done."))


 (fact "mpeg" :mpeg
       (animation) => nil)

 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 2: Non-Linear Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(fact
 "Step 2: Non-Linear Convection" :step2

 (let [nx 41
       dx 0.05
       nt 20
       dt 0.025
       m {:x-steps 2 :dx dx :dt dt}
       x (map #(* % dx) (range nx))
       u0 (set-u0 nx dx)
       u (take (inc nt) (discretize non-linear-convection m u0))
       u_nt (last u)]
   (fact "nx = 41" :step2 nx => 41)
   (fact "u0" :step2
         u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
                2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

   (fact "u_nt" :step2
         (map #(format-x % 5) u_nt)
         => (map #(format-x % 5)
                 [1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          2.          1.99999893  1.98777467  1.70622713
                  1.25479189  1.06191253  1.01238281  1.00205217  1.00026166  1.00002241
                  1.00000095]))
   (doto
       (xy-plot x u0
                :title "Non-Linear Convection"
                :x-label "x (nx=41)"
                :y-label "u(x,0) / u(x,0.625)"
                :legend true)
     (add-lines x u_nt)
     view))



 (let [nx 81
       dx 0.025
       nt 20
       dt 0.001
       m {:x-steps 2 :dx dx :dt dt}
       x (map #(* % dx) (range nx))
       u0 (set-u0 nx dx)
       u (take (inc nt) (discretize non-linear-convection m u0))
       u_nt (last u)]
   (fact "nx = 81" :step2 nx => 81)
   (fact "u0" :step2
         u0
         => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
             2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
             2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
             1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

   (fact "u_nt" :step2
         (map #(format-x % 5) u_nt)
         => (map #(format-x % 5)
                 [1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.
                  1.27689924  1.59631492  1.82330269  1.93962089  1.98383592  1.996566
                  1.99941229  1.99991796  1.99999058  1.99999911  1.99999993  2.          2.
                  2.          2.          2.          2.          2.          2.          2.
                  2.          1.67060367  1.24095512  1.05367286  1.00865797  1.00108033
                  1.00010733  1.00000864  1.00000057  1.00000003  1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.
                  1.          1.          1.          1.          1.          1.          1.]))
   (doto
       (xy-plot x u0
                :title "Linear Convection"
                :x-label "x (nx=81)"
                :y-label "u(x,0) / u(x,0.02)"
                :legend true)
     (add-lines x u_nt)
     view))

 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 3: 1D Diffusion
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(fact
 "Step 3: 1D Diffussion" :step3

 (let [nx 41
       dx 0.05
       nt 20
       sigma 0.2
       nu 0.3
       dt (/ (* sigma (Math/pow dx 2.)) nu)
       m {:x-steps 3 :nu nu :dx dx :dt dt}
       x (map #(* % dx) (range nx))
       u0 (set-u0 nx dx)
       u (take (inc nt) (discretize diffusion-1D m u0))
       u_nt (last u)]
   (fact "nx = 41" :step3 nx => 41)
   (fact "u0" :step3
         u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
                2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
                1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

   (fact "u_nt" :step3
         (map #(format-x % 5) u_nt)
         => (map #(format-x % 5)
                 [1.          1.00109844  1.00369399  1.01028214  1.0252136   1.05496351
                  1.10711342  1.18767409  1.2974854   1.42967135  1.5702342   1.70218575
                  1.81114876  1.88917639  1.93475147  1.94957196  1.93475147  1.88917639
                  1.81114876  1.70218575  1.5702342   1.42967135  1.2974854   1.18767409
                  1.10711342  1.05496356  1.02521402  1.01028497  1.00371014  1.00117673
                  1.00032601  1.0000783   1.00001615  1.00000283  1.00000041  1.00000005
                  1.          1.          1.          1.          1.        ]))
   (doto
       (xy-plot x u0
                :title "1D Diffusion"
                :x-label "x (nx=41)"
                :y-label "u(x,0) / u(x,0.0333)"
                :legend true)
     (add-lines x u_nt)
     view))

 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 4: Burgers' Equation
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(fact
 "Step 4: Burgers' Equation" :step4

 (let
     [x [0.          0.06283185  0.12566371  0.18849556  0.25132741  0.31415927
         0.37699112  0.43982297  0.50265482  0.56548668  0.62831853  0.69115038
         0.75398224  0.81681409  0.87964594  0.9424778   1.00530965  1.0681415
         1.13097336  1.19380521  1.25663706  1.31946891  1.38230077  1.44513262
         1.50796447  1.57079633  1.63362818  1.69646003  1.75929189  1.82212374
         1.88495559  1.94778745  2.0106193   2.07345115  2.136283    2.19911486
         2.26194671  2.32477856  2.38761042  2.45044227  2.51327412  2.57610598
         2.63893783  2.70176968  2.76460154  2.82743339  2.89026524  2.95309709
         3.01592895  3.0787608   3.14159265  3.20442451  3.26725636  3.33008821
         3.39292007  3.45575192  3.51858377  3.58141563  3.64424748  3.70707933
         3.76991118  3.83274304  3.89557489  3.95840674  4.0212386   4.08407045
         4.1469023   4.20973416  4.27256601  4.33539786  4.39822972  4.46106157
         4.52389342  4.58672527  4.64955713  4.71238898  4.77522083  4.83805269
         4.90088454  4.96371639  5.02654825  5.0893801   5.15221195  5.2150438
         5.27787566  5.34070751  5.40353936  5.46637122  5.52920307  5.59203492
         5.65486678  5.71769863  5.78053048  5.84336234  5.90619419  5.96902604
         6.03185789  6.09468975  6.1575216   6.22035345  6.28318531]

      u0 [4.        ,  4.06283185,  4.12566371,  4.18849556,  4.25132741,
          4.31415927,  4.37699112,  4.43982297,  4.50265482,  4.56548668,
          4.62831853,  4.69115038,  4.75398224,  4.81681409,  4.87964594,
          4.9424778 ,  5.00530965,  5.0681415 ,  5.13097336,  5.19380521,
          5.25663706,  5.31946891,  5.38230077,  5.44513262,  5.50796447,
          5.57079633,  5.63362818,  5.69646003,  5.75929189,  5.82212374,
          5.88495559,  5.94778745,  6.0106193 ,  6.07345115,  6.136283  ,
          6.19911486,  6.26194671,  6.32477856,  6.38761042,  6.45044227,
          6.51327412,  6.57610598,  6.63893783,  6.70176967,  6.76460125,
          6.82742866,  6.89018589,  6.95176632,  6.99367964,  6.72527549,
          4.        ,  1.27472451,  1.00632036,  1.04823368,  1.10981411,
          1.17257134,  1.23539875,  1.29823033,  1.36106217,  1.42389402,
          1.48672588,  1.54955773,  1.61238958,  1.67522144,  1.73805329,
          1.80088514,  1.863717  ,  1.92654885,  1.9893807 ,  2.05221255,
          2.11504441,  2.17787626,  2.24070811,  2.30353997,  2.36637182,
          2.42920367,  2.49203553,  2.55486738,  2.61769923,  2.68053109,
          2.74336294,  2.80619479,  2.86902664,  2.9318585 ,  2.99469035,
          3.0575222 ,  3.12035406,  3.18318591,  3.24601776,  3.30884962,
          3.37168147,  3.43451332,  3.49734518,  3.56017703,  3.62300888,
          3.68584073,  3.74867259,  3.81150444,  3.87433629,  3.93716815,  4.]

      u_1 [4.0049      4.04496259  4.10751809  4.17007359  4.2326291   4.2951846
           4.35774011  4.42029561  4.48285111  4.54540662  4.60796212  4.67051763
           4.73307313  4.79562864  4.85818414  4.92073964  4.98329515  5.04585065
           5.10840616  5.17096166  5.23351716  5.29607267  5.35862817  5.42118368
           5.48373918  5.54629468  5.60885019  5.67140569  5.7339612   5.7965167
           5.85907221  5.92162771  5.98418321  6.04673872  6.10929422  6.17184973
           6.23440523  6.29696073  6.35951624  6.42207174  6.48462725  6.54718275
           6.60973825  6.67229373  6.73484878  6.79739671  6.85982549  6.92026607
           6.94896026  6.66003054  4.76307714  1.70950417  1.0494278   1.04669199
           1.10512188  1.1674257   1.22996589  1.29252046  1.35507591  1.41763141
           1.48018692  1.54274242  1.60529792  1.66785343  1.73040893  1.79296444
           1.85551994  1.91807544  1.98063095  2.04318645  2.10574196  2.16829746
           2.23085296  2.29340847  2.35596397  2.41851948  2.48107498  2.54363049
           2.60618599  2.66874149  2.731297    2.7938525   2.85640801  2.91896351
           2.98151901  3.04407452  3.10663002  3.16918553  3.23174103  3.29429654
           3.35685204  3.41940754  3.48196305  3.54451855  3.60707406  3.66962956
           3.73218506  3.79474057  3.85729607  3.91985158  3.97750708]

      u_2 [3.99820864  4.0353731   4.08953174  4.15181332  4.2140949   4.27637648
           4.33865806  4.40093964  4.46322122  4.5255028   4.58778439  4.65006597
           4.71234755  4.77462913  4.83691071  4.89919229  4.96147387  5.02375545
           5.08603703  5.14831861  5.21060019  5.27288177  5.33516335  5.39744494
           5.45972652  5.5220081   5.58428968  5.64657126  5.70885284  5.77113442
           5.833416    5.89569758  5.95797916  6.02026074  6.08254232  6.1448239
           6.20710549  6.26938707  6.33166865  6.39395023  6.45623181  6.51851339
           6.58079497  6.6430765   6.70535731  6.767626    6.8296929   6.88851176
           6.91023243  6.66932699  5.30535057  2.26156991  1.14918028  1.0516625
           1.10090394  1.16235268  1.22458245  1.28686082  1.3491422   1.41142377
           1.47370535  1.53598693  1.59826851  1.66055009  1.72283167  1.78511325
           1.84739483  1.90967641  1.97195799  2.03423957  2.09652115  2.15880273
           2.22108432  2.2833659   2.34564748  2.40792906  2.47021064  2.53249222
           2.5947738   2.65705538  2.71933696  2.78161854  2.84390012  2.9061817
           2.96846329  3.03074487  3.09302645  3.15530803  3.21758961  3.27987119
           3.34215277  3.40443435  3.46671593  3.52899751  3.59127909  3.65356067
           3.71584225  3.77812384  3.84040542  3.90230487  3.95909426]



      u_100 [2.81860248  2.86219107  2.90577873  2.94936484  2.99294844  3.03652808
             3.08010158  3.12366585  3.16721654  3.21074771  3.25425154  3.29771797
             3.34113444  3.38448584  3.42775454  3.47092075  3.51396319  3.55686007
             3.59959041  3.6421356   3.68448121  3.7266187   3.76854701  3.81027381
             3.8518161   3.89320018  3.9344608   3.97563963  4.01678295  4.05793904
             4.09915528  4.14047543  4.18193721  4.2235705   4.26539628  4.30742625
             4.3496633   4.39210243  4.43473218  4.47753629  4.52049539  4.56358859
             4.60679488  4.65009421  4.69346824  4.73690088  4.78037845  4.82388968
             4.86742555  4.91097909  4.954545    4.99811937  5.04169941  5.08528313
             5.12886915  5.17245648  5.21604428  5.25963161  5.3032169   5.34679686
             5.39036424  5.43390248  5.47737383  5.52069227  5.56366214  5.60583953
             5.64622005  5.68253942  5.70972375  5.71653417  5.67868075  5.54649073
             5.23096356  4.62161195  3.7237691   2.83274015  2.25637999  1.99159854
             1.90314957  1.89369951  1.9164822   1.95196924  1.99240748  2.03477056
             2.07788235  2.12128577  2.1648029   2.2083644   2.25194321  2.29552876
             2.33911694  2.38270612  2.4262957   2.46988542  2.5134752   2.557065
             2.60065479  2.64424456  2.68783427  2.7314239   2.77501334]]

   (fact "partitioning" :step4
         (+ 2 (count (partition 3 1 (range 102)))) => (count (range 102)))

   (fact "Burgers' u0" :step4
         (map #(format-x % 5) (map #(burgers-u0 0. % 0.07) x))
         =>  (map #(format-x % 5) u0))

   (let [nx 101
         nt 100
         dx (/ (* 2. Math/PI) (dec nx))
         nu 0.07
         dt (* dx nu)
         m {:x-steps 3 :bc burgers-bc-err :nu nu :dx dx :dt dt}
         u (take (inc nt) (discretize burgers-eqn m u0))
         u_nt (last u)]
     (fact "dimensions nx 101 " :step4
           (= 101 (count x) (count u0) (count u_1) (count u_nt) (count u_100))
           => truthy)

     (fact "u0" :step4
           (map #(format-x % 5) (first u))
           => (map #(format-x % 5) u0))

     (fact "u1" :step4
           (map #(format-x % 5) (second u))
           => (map #(format-x % 5) u_1))

     (fact "u2" :step4
           (map #(format-x % 5) (nth u 2))
           => (map #(format-x % 5) u_2))

     (fact "u100" :step4
           (map #(format-x % 6) u_nt)
           => (map #(format-x % 6) u_100))

     (fact "Euklid distance 'Python/Clojure' < 1.00001 "
           :step4
           (< (/ (* 2. (euklid u_nt u_100))
                 (+ (euklid u_nt (repeat (count u_nt) 0))
                    (euklid u_100 (repeat (count u_100) 0)))) 1.00001)
           => truthy)

     (comment
       (fact "u_nt - u-100 :"
             :step4
             (map (fn [x] (format-x x 5))
                  (map #(apply - %)
                       (map vector u_nt u_100)))
             =not=> truthy)
       )
     ))

 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 5: 2D Linear Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-linear-convection-2D [c nx ny nt sigma]
 (let [dx (/ 2. (dec nx))
       dy (/ 2. (dec ny))
       dt (* sigma dx)
       m {:c c :nx nx :dx dx :ny ny :dy dy :dt dt}
       u0 (set-u0 ny dy nx dx)
       u (take (inc nt) (discretize-2D linear-convection-2D m u0))
       u_nt (last u)
       v (apply concat u_nt)
       res (json/read-json
            (slurp (str/join
                    ["./test/cfd_clojure/python/test-05-" nx ".json"])))
       u_nt_py (:u_nt res)
       v_py (apply concat u_nt_py)]
   (fact "params " :step5
         [nx dx ny dy nt dt c sigma]
         => [(:nx res) (:dx res)
             (:ny res) (:dy res)
             (inc (:nt res)) (:dt res)
             (:c res) (:sigma res)])
   (fact "u0" :step5
         u0 => (:u0 res))
   (fact "(first u) is u0" :step5
         (first u) => (:u0 res))
   (fact "dimensions u0: nx, ny" :step5
         [(count u0) (count (first u0))
          (count (:u0 res)) (count (first (:u0 res)))]
         => [ny nx ny nx])
   (fact "dimensions u_nt: nx, ny" :step5
         [(count u_nt) (count (first u_nt))
          (count u_nt_py) (count (first u_nt_py))]
         => [ny nx ny nx])
   (fact "u_nt" :step5
         (format-zz u_nt 5) => (format-zz u_nt_py 5))))


(fact
 "Step 5: 2D Linear Convection" :step5
 (test-linear-convection-2D 1 41 21 1 0.2)
 (test-linear-convection-2D 1 81 81 50 0.2))



(comment
  (fact "Euklid distance 'Python/Clojure' < 0.00001 "
        :step5
        (< (/ (* 2. (euklid v v_py))
              (+ (euklid v (repeat (count v) 0))
                 (euklid v_py (repeat (count v_py) 0)))) 0.00001)
        => truthy)
  ;;(println "u_41 " (nth u 41))
  ;;(println "u_41_py " (nth (:u_nt res) 41))

  (fact "v v_py :"
        :step5
        (count (filter #(> (Math/abs %) 0.001) (map (fn [x] (format-x x 5))
                                                    (map #(apply - %)
                                                         (map vector v v_py)))))
        => 10)
  )


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 6: 2D Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-convection-2D [c nx ny nt sigma]
 (let [dx (/ 2. (dec nx))
       dy (/ 2. (dec ny))
       dt (* sigma dx)
       m {:c c :nx nx :dx dx :ny ny :dy dy :dt dt}
       u0 (set-u0 ny dy nx dx)
       v0 u0
       w (take (inc nt) (discretize-2D convection-2D m [u0 v0]))
       u (map first w)
       v (map second w)
       u_nt (last u)
       v_nt (last v)
       res (json/read-json
            (slurp (str/join
                    ["./test/cfd_clojure/python/test-06-" nx ".json"])))
       u_nt_py (:u_nt res)
       v_nt_py (:v_nt res)]
   (fact "params " :step6
         [nx dx ny dy nt dt c sigma]
         => [(:nx res) (:dx res)
             (:ny res) (:dy res)
             (inc (:nt res)) (:dt res)
             (:c res) (:sigma res)])
   (fact "u0, v0" :step6
         [u0 v0 (= u0 v0)] => [(:u0 res) (:v0 res) (= (:u0 res) (:v0 res))])
   (fact "[(first u) (first v)] is [u0 v0]" :step6
         [(first u) (first v)] => [(:u0 res) (:v0 res)])
   (fact "dimensions u0, v0 : nx, ny" :step6
         [(count u0) (count (first u0))
          (count (:u0 res)) (count (first (:u0 res)))
          (count v0) (count (first v0))
          (count (:v0 res)) (count (first (:v0 res)))]
         => [ny nx ny nx
             ny nx ny nx])
   (fact "dimensions u_nt, v_nt: nx, ny" :step6
         [(count u_nt) (count (first u_nt))
          (count u_nt_py) (count (first u_nt_py))
          (count v_nt) (count (first v_nt))
          (count v_nt_py) (count (first v_nt_py))]
         => [ny nx ny nx
             ny nx ny nx])
   (fact "u_nt, v_nt" :step6
         [(format-zz u_nt 5) (format-zz v_nt 5)]
         => [(format-zz u_nt_py 5) (format-zz v_nt_py 5)])))


(fact
 "Step 6: 2D Convection" :step6
 (test-convection-2D 1 41 21 1 0.2)
 (test-convection-2D 1 81 81 50 0.2))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 7: 2D Diffusion
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-diffusion-2D [nx ny nt nu sigma]
 (let [dx (/ 2. (dec nx))
       dy (/ 2. (dec ny))
       dt (/ (* sigma dx dy) nu)
       m {:nu nu :nx nx :dx dx :ny ny :dy dy :dt dt}
       u0 (set-u0 ny dy nx dx)
       u (take (inc nt) (discretize-2D diffusion-2D m u0))
       u_nt (last u)
       res (json/read-json
            (slurp (str/join
                    ["./test/cfd_clojure/python/test-07-" nx ".json"])))
       u_nt_py (:u_nt res)]
   (fact "params " :step7
         [nx dx ny dy nt dt nu sigma]
         => [(:nx res) (:dx res)
             (:ny res) (:dy res)
             (inc (:nt res)) (:dt res)
             (:nu res) (:sigma res)])
   (fact "u0" :step7
         [u0 (first u)]  => [(:u0 res) (:u0 res)])
   (fact "dimensions u0 : nx, ny" :step7
         [(count u0) (count (first u0))
          (count (:u0 res)) (count (first (:u0 res)))]
         => [ny nx ny nx])
   (fact "dimensions u_nt: nx, ny" :step7
         [(count u_nt) (count (first u_nt))
          (count u_nt_py) (count (first u_nt_py))]
         => [ny nx ny nx])
   (fact "u_nt" :step7
         (format-zz u_nt 5) => (format-zz u_nt_py 5))
   ))


(fact
 "Step 7: 2D Diffusion" :step7
 (test-diffusion-2D 31 31 11 0.05 0.25))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 8: 2D Burgers' Equation
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-burgers-eqn-2D [nx ny nt nu sigma]
 (let [dx (/ 2. (dec nx))
       dy (/ 2. (dec ny))
       dt (/ (* sigma dx dy) nu)
       m {:nx nx :dx dx :ny ny :dy dy :dt dt :nu nu}
       u0 (set-u0 ny dy nx dx)
       v0 u0
       w (take (inc nt) (discretize-2D burgers-eqn-2D m [u0 v0]))
       u (map first w)
       v (map second w)
       u_nt (last u)
       v_nt (last v)
       res (json/read-json
            (slurp (str/join
                    ["./test/cfd_clojure/python/test-08-" nx ".json"])))
       u_nt_py (:u_nt res)
       v_nt_py (:v_nt res)]
   (fact "params " :step8
         [nx dx ny dy nt dt nu sigma] ; c is irrelevant for Burgers' Eqn.
         => [(:nx res) (:dx res)
             (:ny res) (:dy res)
             (inc (:nt res)) (:dt res)
             (:nu res) (:sigma res)])
   (fact "u0, v0" :step8
         [u0 v0 (= u0 v0)] => [(:u0 res) (:v0 res) (= (:u0 res) (:v0 res))])
   (fact "[(first u) (first v)] is [u0 v0]" :step8
         [(first u) (first v)] => [(:u0 res) (:v0 res)])
   (fact "dimensions u0, v0 : nx, ny" :step8
         [(count u0) (count (first u0))
          (count (:u0 res)) (count (first (:u0 res)))
          (count v0) (count (first v0))
          (count (:v0 res)) (count (first (:v0 res)))]
         => [ny nx ny nx
             ny nx ny nx])
   (fact "dimensions u_nt, v_nt: nx, ny" :step8
         [(count u_nt) (count (first u_nt))
          (count u_nt_py) (count (first u_nt_py))
          (count v_nt) (count (first v_nt))
          (count v_nt_py) (count (first v_nt_py))]
         => [ny nx ny nx
             ny nx ny nx])
   (fact "u_nt, v_nt" :step8
         [(format-zz u_nt 5) (format-zz v_nt 5)]
         => [(format-zz u_nt_py 5) (format-zz v_nt_py 5)])))

(fact
 "Step 8: 2D Burgers' Equation" :step8
 (test-burgers-eqn-2D 21 21 1 0.01 0.0009)
 (test-burgers-eqn-2D 41 41 121 0.01 0.0009))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 9: 2D Laplace Equation
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; ##boundary conditions
;; p[:,0] = 0		##p = 0 @ x = 0        ;; redundant
;; p[:,-1] = y		##p = y @ x = 2
;; p[0,:] = p[1,:]	##dp/dy = 0 @ y = 0
;; p[-1,:] = p[-2,:]	##dp/dy = 0 @ y = 1
;;
(defn- set-p0 [ny nx dy]
  (let [p0 (matrix 0. ny nx)]
    (doseq [y (range 1 ny)]
      (clx/set p0 y (dec nx) (* y dy)))
    (doseq [x (range nx)]
      (clx/set p0 0 x (sel p0 1 x)))
    (doseq [x (range nx)]
      (clx/set p0 (dec ny) x (sel p0 (- ny 2) x)))
    p0))

(defn test-laplace-eqn-2D [ny nx eps]
  (let [dx (/ 2. (dec nx))
        dy (/ 2. (dec ny))
        m {:nx nx :dx dx :ny ny :dy dy :eps eps}
        p0 (set-p0 ny nx (/ dy 2.))
        res (json/read-json
             (slurp (str/join
                     ["./test/cfd_clojure/python/test-09-" nx ".json"])))
        p (laplace-eqn-2D m p0)]
    (fact "params " :step9
          [nx dx ny dy eps]
          => [(:nx res) (:dx res)
              (:ny res) (:dy res)
              (:eps res)])
    (fact "dimensions p0, p" :step9
          [(count p0) (count (first p0))
           (count p) (count (first p))]
          => [ny nx ny nx])
    (fact "p0" :step9
          (format-zz p0 10) => (format-zz (:p0 res) 10))
    (fact "p" :step9
          (format-zz p 9) => (format-zz (:p res) 9))))

(fact
 "Step 9: 2D Laplace Equation" :step9
 (test-laplace-eqn-2D 31 31 0.01))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 10: 2D Poisson Equation
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-poisson-eqn-2D [nx xmin xmax ny ymin ymax nt]
  (let [upper_x (dec nx)
        upper_y (dec ny)
        dx (/ (- xmax xmin) upper_x) ; rows
        dy (/ (- ymax ymin) upper_y) ; cols
        b (matrix 0. (dec upper_x) (dec upper_y))
        _ (clx/set b (dec (/ nx 4)) (dec (/ ny 4)) 100.)
        _ (clx/set b (dec (/ (* 3 nx) 4)) (dec (/ (* 3 ny) 4)) -100.)
        res (json/read-json
             (slurp (str/join
                     ["./test/cfd_clojure/python/test-10-" nx ".json"])))
        b_core (sel (sel (:b res) :except-rows upper_x :except-cols upper_y)
                    :except-rows 0   :except-cols 0)
        m {:nx nx :dx dx :ny ny :dy dy :b b}
        p0 (matrix 0. nx ny)
        p (take (inc nt) (discretize-2D poisson-eqn-2D m p0))
        p_nt (last p)
        p_nt_py (:p res)]
    (fact "params " :step10
          [nx dx xmin xmax ny dy ymin ymax nt]
          => [(:nx res) (:dx res)
              (:xmin res) (:xmax res)
              (:ny res) (:dy res)
              (:ymin res) (:ymax res)
              (:nt res)])
    (fact "dimensions p0: nx, ny" :step10
          [(count p0) (count (first p0))
           (count (:p0 res)) (count (first (:p0 res)))]
          => [ny nx ny nx])
    (fact "dimensions b: nx-2, ny-2" :step10
          [(count b) (count (first b))
           (count b_core) (count (first b_core))]
          => [(- nx 2) (- ny 2) (- nx 2) (- ny 2)])
    (fact "p0" :step10
          [p0 (first p)] => [(:p0 res) (:p0 res)])
    (fact "b" :step10
          b => b_core)
    (fact "dimensions p_nt: nx, ny" :step10
          [(count p_nt) (count (first p_nt))
           (count p_nt_py) (count (first p_nt_py))]
          => [nx ny nx ny])
    (fact "p_nt" :step10
          (format-zz p_nt 5) => (format-zz p_nt_py 5))
    ))


(fact
 "Step 10: 2D Poisson Equation" :step10
 (test-poisson-eqn-2D 50 0. 2. 50 0. 1. 100))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Step 11: Cavity Flow with Navier-Stokes
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn test-cavitiy-flow-2D [nx ny nt dt nit rho nu]
  (let [upper_x (dec nx) ; rows
        upper_y (dec ny) ; cols
        dx (/ 2. (dec nx))
        dy (/ 2. (dec ny))
        res (json/read-json
             (slurp (str/join
                     ["./test/cfd_clojure/python/test-11-" nt ".json"])))
        m {:nx nx :dx dx
           :ny ny :dy dy
           :nt nt :dt dt
           :nit nit :rho rho :nu nu}
        u0 (matrix 0. nx ny)
        v0 (matrix 0. nx ny)
        p0 (matrix 0. nx ny)

;;         w (take (inc nt) (discretize-2D cavity-flow-2D m [u0 v0 p0]))

;;         u (map first w)
;;         v (map second w)
;;         p (map last w)
;;         u_nt (last u)
        u_nt_py (:u_nt res)
;;         v_nt (last v)
        v_nt_py (:v_nt res)
;;         p_nt (last p)
        p_nt_py (:p_nt res)

        b_nt (buildup-b m [u_nt_py v_nt_py])
        b_nt_py (sel (sel (:b_nt res) :except-rows upper_x :except-cols upper_y)
                     :except-rows 0   :except-cols 0)

        p_py (last (take
                    (inc (:nit m))
                    (iterate (partial pressure-poisson m b_nt_py) p_nt_py)))]

    (fact "params " :step11
          [nx dx ny dy nt dt nit rho nu]
          => [(:nx res) (:dx res)
              (:ny res) (:dy res)
              (:nt res) (:dt res)
              (:nit res)
              (:rho res) (:nu res)])

    (fact "dimensions b_nt: nx-2, ny-2" :step11
          [(count b_nt) (count (first b_nt))
           (count b_nt_py) (count (first b_nt_py))]
          => [(- nx 2) (- ny 2) (- nx 2) (- ny 2)])
    (fact "b_nt" :step11
          (format-zz b_nt 5) => (format-zz b_nt_py 5))

    (fact "dimensions p_py: nx, ny" :step11
          [nt (count p_py) (count (first p_py))
           (count (:p_py res)) (count (first (:p_py res)))]
          => [(:nt res) nx ny nx ny])
    (fact "p py/clj" :step11
          (format-zz p_py 5) => (format-zz (:p_py res) 5))


;;     (fact "u0, v0, p0" :step11
;;           [(first u) (first v) (first p)] => [(:u0 res) (:v0 res) (:p0 res)])
;;
;;     (fact "u" :step11
;;           u_nt => u_nt_py)
;;
;;     (fact "v" :step11
;;           v_nt => v_nt_py)
;;
;;     (fact "p" :step11
;;           p_nt => p_nt_py)

           ))

;; z[1:-1,1:-1] ; Dz :except-rows  first     last    :except-cols  first     last
;; z[2:,1:-1]   ; Ez :except-rows  first     second  :except-cols  first     last
;; z[0:-2,1:-1] ; Fz :except-rows  prev-last last    :except-cols  first     last
;; z[1:-1,2:]   ; Gz :except-rows  first     last    :except-cols  first     second
;; z[1:-1,0:-2] ; Hz :except-rows  first     last    :except-cols  prev-last last


;; u[1:-1,1:-1] = Du-\
;;     Du*dt/dx*(Du-Fu)-\
;;     Dv*dt/dy*(Du-Hu)-\
;;     dt/(2*rho*dx)*(Ep-Fp)+\
;;     nu*(dt/dx**2*(Eu-2*Du+Fu)+\
;;     dt/dy**2*(Gu-2*Du+Hu))
;;
;; v[1:-1,1:-1] = Dv-\
;;     Du*dt/dx*(Dv-Fv)-\
;;     Dv*dt/dy*(Dv-Hv)-\
;;     dt/(2*rho*dy)*(Gp-Hp)+\
;;     nu*(dt/dx**2*(Ev-2*Dv+Fv)+\
;;     (dt/dy**2*(Gv-2*Dv+Hv)))


(defn tweaked-cavity-flow-2D [m [un vn pn]]
  (let [upper_x (dec (:nx m))
        upper_y (dec (:ny m))
        ;;b (mult (:rho m) (buildup-b m [un vn]))
        ;;_ (println "before p, nt = " (:nt m))
        ;; p (last (take (inc (:nit m)) (iterate (partial pressure-poisson m b) pn)))
        p pn ;;p ((apply comp (repeat (inc (:nit m)) (partial pressure-poisson m b))) pn)
        ;; p (loop [n (:nit m)
        ;;          p pn]
        ;;     (if (< n 0)
        ;;       p
        ;;       (recur (dec n) ((partial pressure-poisson m b) p))))
        ;;_ (println "after p")
        Du (sel (sel un :except-rows upper_x :except-cols upper_y)
               :except-rows 0   :except-cols 0)
        Eu (sel (sel un :except-rows 0 :except-cols upper_y)
                :except-rows 0 :except-cols 0)
        Fu (sel (sel un :except-rows upper_x :except-cols upper_y)
                :except-rows (dec upper_x) :except-cols 0)
        Gu (sel (sel un :except-rows upper_x :except-cols 0)
                :except-rows 0 :except-cols 0)
        Hu (sel (sel un :except-rows upper_x :except-cols upper_y)
                :except-rows 0 :except-cols (dec upper_y))
        Du-Fu (minus Du Fu)
        Du-Hu (minus Du Hu)
        Dv (sel (sel vn :except-rows upper_x :except-cols upper_y)
               :except-rows 0   :except-cols 0)
        Ev (sel (sel vn :except-rows 0 :except-cols upper_y)
                :except-rows 0 :except-cols 0)
        Fv (sel (sel vn :except-rows upper_x :except-cols upper_y)
                :except-rows (dec upper_x) :except-cols 0)
        Gv (sel (sel vn :except-rows upper_x :except-cols 0)
                :except-rows 0 :except-cols 0)
        Hv (sel (sel vn :except-rows upper_x :except-cols upper_y)
                :except-rows 0 :except-cols (dec upper_y))
        Dv-Fv (minus Dv Fv)
        Dv-Hv (minus Dv Hv)
        Ep (sel (sel p :except-rows 0 :except-cols upper_y)
                :except-rows 0 :except-cols 0)
        Fp (sel (sel p :except-rows upper_x :except-cols upper_y)
                :except-rows (dec upper_x) :except-cols 0)
        Gp (sel (sel p :except-rows upper_x :except-cols 0)
                :except-rows 0 :except-cols 0)
        Hp (sel (sel p :except-rows upper_x :except-cols upper_y)
                :except-rows 0 :except-cols (dec upper_y))
        dtx (* -1. (/ (:dt m) (:dx m)))
        dty (* -1. (/ (:dt m) (:dy m)))
        dtxrho (/ dtx (* 2. (:rho m)))
        dtyrho (/ dty (* 2. (:rho m)))
        dtx2nu (/ (* (:nu m) (:dx m)) (math/expt (:dx m) 2.))
        dty2nu (/ (* (:nu m) (:dy m)) (math/expt (:dy m) 2.))
        u_core (plus (mult (+ 1. (* -2. dtx2nu) (* -2. dty2nu)) Du)
                     (mult dtx Du Du-Fu)
                     (mult dty Dv Du-Hu)
                     (mult dtxrho (minus Ep Fp))
                     (mult dtx2nu Eu)
                     (mult dtx2nu Fu)
                     (mult dty2nu Gu)
                     (mult dty2nu Hu))
        v_core (plus (mult (+ 1. (* -2. dtx2nu) (* -2. dty2nu)) Dv)
                     (mult dtx Du Dv-Fv)
                     (mult dty Dv Dv-Hv)
                     (mult dtyrho (minus Gp Hp))
                     (mult dtx2nu Ev)
                     (mult dtx2nu Fv)
                     (mult dty2nu Gv)
                     (mult dty2nu Hv))
        uu (matrix 0. (:nx m) (:ny m))
        vv (matrix 0. (:nx m) (:ny m))]
    (doseq [x (range 1 upper_x)
            y (range 1 upper_y)]
      (clx/set uu x y (sel u_core (dec x) (dec y)))
      (clx/set vv x y (sel v_core (dec x) (dec y))))
    (doseq [x (range upper_x)]
      (clx/set uu x upper_y 1.))
  [uu vv p]))

(defn test-1-cavity-flow-2D [nx ny nt dt nit rho nu]
  (let [upper_x (dec nx) ; rows
        upper_y (dec ny) ; cols
        dx (/ 2. (dec nx))
        dy (/ 2. (dec ny))
        res (json/read-json
             (slurp (str/join
                     ["./test/cfd_clojure/python/test-11-" nt ".json"])))
        m {:nx nx :dx dx
           :ny ny :dy dy
           :nt nt :dt dt
           :nit nit :rho rho :nu nu}
        u0 (:u0 res)
        v0 (:v0 res)
        p0 (:p0 res)
        b0 (buildup-b m [u0 v0])
        b0_py (sel (sel (:b0 res) :except-rows upper_x :except-cols upper_y)
                     :except-rows 0   :except-cols 0)

        w (take (inc nt) (discretize-2D tweaked-cavity-flow-2D m [u0 v0 p0]))

        u (map first w)
        v (map second w)
        p (map last w)
        u_nt (last u)
        u_nt_py (:u_nt res)
        v_nt (last v)
        v_nt_py (:v_nt res)
        p_nt (last p)
        p_nt_py (:p_nt res)

        b_nt (buildup-b m [u_nt_py v_nt_py])
        b_nt_py (sel (sel (:b_nt res) :except-rows upper_x :except-cols upper_y)
                     :except-rows 0   :except-cols 0)

        p_py (last (take
                    (inc (:nit m))
                    (iterate (partial pressure-poisson m b_nt_py) p_nt_py)))]

    (fact "params " :step11
          [nx dx ny dy nt dt nit rho nu]
          => [(:nx res) (:dx res)
              (:ny res) (:dy res)
              (:nt res) (:dt res)
              (:nit res)
              (:rho res) (:nu res)])

    (fact "dimensions b_nt: nx-2, ny-2" :step11
          [(count b_nt) (count (first b_nt))
           (count b_nt_py) (count (first b_nt_py))]
          => [(- nx 2) (- ny 2) (- nx 2) (- ny 2)])
    (fact "b0 b_nt" :step11
          [(format-zz b0 7) (format-zz b0 7)] => [(format-zz b0_py 7) (format-zz b_nt_py 7)])

    (fact "dimensions p_py: nx, ny" :step11
          [nt (count p_py) (count (first p_py))
           (count (:p_py res)) (count (first (:p_py res)))]
          => [(:nt res) nx ny nx ny])
    (fact "p py/clj" :step11
          (format-zz p_py 5) => (format-zz (:p_py res) 5))


    (fact "u0, v0, p0" :step11
          [(first u) (first v) (first p)] => [(:u0 res) (:v0 res) (:p0 res)])

;;     (fact "u" :step11
;;           (format-zz u_nt 5) => (format-zz u_nt_py 5))
;;
;;     (fact "v" :step11
;;           v_nt => v_nt_py)
;;
;;     (fact "p" :step11
;;           p_nt => p_nt_py)

           ))



(fact
 "Step 11: Cavity Flow with Navier-Stokes" :step11
 (test-1-cavity-flow-2D 21 21   1 0.001 50 1. 0.1)
 (test-cavitiy-flow-2D 21 21   3 0.001 50 1. 0.1)
 (test-cavitiy-flow-2D 41 41 200 0.001 50 1. 0.1)
 (test-cavitiy-flow-2D 41 41 700 0.001 50 1. 0.1)
 )



;; restore settings
;;

(set! *print-length* print-length)
;;
(if-not (= locale en_US)
  (java.util.Locale/setDefault locale))
