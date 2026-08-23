import { dirname } from "path";
import { fileURLToPath } from "url";
import { FlatCompat } from "@eslint/eslintrc";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const compat = new FlatCompat({
  baseDirectory: __dirname,
});

const eslintConfig = [
  // `next lint` applied these implicitly. The lint script now calls `eslint .`
  // directly - `next lint` is deprecated in Next 15.5, and it was passing
  // eslintrc-era options to the flat-config engine, so it could not start at
  // all - which means the ignores have to be stated here.
  {
    ignores: [
      ".next/**",
      "out/**",
      "build/**",
      "dist/**",
      "next-env.d.ts",
      "electron/dist/**",
      "electron/backend/**",
      "electron/standalone/**",
    ],
  },
  ...compat.extends("next/core-web-vitals", "next/typescript"),
  {
    // Electron main/preload and the build helper are plain CommonJS run by
    // Node, not bundled modules - require() is the correct form there.
    files: ["electron/**/*.js", "scripts/**/*.js"],
    rules: {
      "@typescript-eslint/no-require-imports": "off",
    },
  },
];

export default eslintConfig;
