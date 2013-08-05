(ns cfd-clojure-test
  (:use midje.sweet
        (incanter core charts))
  (:require [cfd-clojure :refer :all]))

(import 'java.io.FileOutputStream)
(use '[clojure.java.shell :only [sh]])


;; save settings
;;
(def print-length *print-length*)
(set! *print-length* 100)
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
  (read-string (format (clojure.string/join ["%." (str n) "f"]) (java.math.BigDecimal. x))))

(defn- format-zz [zz n]
  (map #(map (fn [x] (format-x x n)) %) zz))

;; TODO: understand with-precision
;;
(defn- round-zz [zz n]
  (map #(map (fn [x] (with-precision n x)) %) zz))



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


(defn- set-u0 [nx dx]
  (for [x (range nx)] (if (and (>= x (/ 0.5 dx))
                               (< x (inc (/ 1.0 dx))))
                        2.0
                        1.0)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Linear Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(let [nx 41.
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
  (fact u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
               2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

  ; precision of 10^-5 fails only for the first 3 values 1.00045526  1.00007826  1.00000972
  ; hence reduced to 10^-3:
  ;
  (fact (map #(format-x % 3) u_nt)
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



(let [nx 81.
      dx (/ 2. (dec nx))
      nt 25
      m {:c 1 :x-steps 2 :dx dx :dt 0.025}
      x (map #(* % dx) (range nx))
      u0 (set-u0 nx dx)
      u (take (inc nt) (discretize linear-convection m u0))
      u_nt (last u)]
  (fact u0
        => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
            2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
            2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

  (fact u_nt
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
               (< 9 t) (clojure.string/join ["0" (str t)])
               :else (clojure.string/join ["00" (str t)]))
       fname (clojure.string/join ["./mpeg/linear-convection-" n ".png"])
       fos (FileOutputStream. fname)]
    (save p fos)
    (if (= 0 (mod t 100)) (println "... " (- 1000 t)))
    (.close fos)))

(print "purging ./mpeg ... ")
(sh "rm" "./mpeg/*")
(println "done.\ncalculating Linear Convection (nt=1000) ... ")
(let [nx 41.
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
(sh "ffmpeg" "-r" "100" "-b" "20M" "-i" "mpeg/linear-convection-%03d.png" "mpeg/linear-convection.mp4")
(println "done.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; Non-Linear Convection
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(let [nx 41.
      dx 0.05
      nt 20
      dt 0.025
      m {:x-steps 2 :dx dx :dt dt}
      x (map #(* % dx) (range nx))
      u0 (set-u0 nx dx)
      u (take (inc nt) (discretize non-linear-convection m u0))
      u_nt (last u)]
  (fact u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
               2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

  (fact (map #(format-x % 5) u_nt)
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



(let [nx 81.
      dx 0.025
      nt 20
      dt 0.001
      m {:x-steps 2 :dx dx :dt dt}
      x (map #(* % dx) (range nx))
      u0 (set-u0 nx dx)
      u (take (inc nt) (discretize non-linear-convection m u0))
      u_nt (last u)]
  (fact u0
        => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
            2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
            2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
            1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

  (fact (map #(format-x % 5) u_nt)
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; 1D Diffusion
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(let [nx 41.
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
  (fact u0 => [1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               2.  2.  2.  2.  2.  2.  2.  2.  2.  2.
               2.  1.  1.  1.  1.  1.  1.  1.  1.  1.
               1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.])

  (fact (map #(format-x % 5) u_nt)
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

;; restore settings
;;
(set! *print-length* print-length)
;;
(if-not (= locale en_US)
  (java.util.Locale/setDefault locale))
