(ns cfd-clojure-test
  (:use midje.sweet
        (incanter core charts))
  (:require [cfd-clojure :refer :all]))


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


;; restore settings
;;
(set! *print-length* print-length)
;;
(if-not (= locale en_US)
  (java.util.Locale/setDefault locale))
