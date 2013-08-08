(defproject cfd-clojure "0.0.1-SNAPSHOT"
  :description "CFD Python port to Clojure"
  :dependencies [[org.clojure/clojure "1.5.1"]

                 [org.clojure/math.numeric-tower "0.0.2"]

                 [org.clojure/core.match "0.2.0-rc5"]

                 [incanter/incanter-core "1.5.2"]
                 ;[incanter/incanter-io "1.5.2"]
                 [incanter/incanter-charts "1.5.2"]
                 ;[incanter/incanter-mongodb "1.5.2"]
                 ;[incanter/incanter-pdf "1.5.2"]
                 ;[incanter/incanter-latex "1.5.2"]
                 ;[incanter/incanter-excel "1.5.2"]
                 ;[incanter/incanter-sql "1.5.2"]
                 ;[incanter/incanter-zoo "1.5.2"]
                 [swingrepl "1.3.0"
                  :exclusions [org.clojure/clojure
                               org.clojure/clojure-contrib]]
                 [jline "0.9.94"]]

  :profiles {:dev {:dependencies [[midje "1.5.1"]]}}
  :repl-options {:port 4555})
