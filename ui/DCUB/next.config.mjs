/** @type {import('next').NextConfig} */
const nextConfig = {
  // output: 'export',
  output: 'standalone',
  trailingSlash: true,
  // Both nets are enforced again - the build now reports "Linting and checking
  // validity of types". While they were suppressed, `npm run lint` could not
  // even start (it invoked the deprecated `next lint`, which handed
  // eslintrc-era options to the flat-config engine) and a TS2339 in
  // lib/electron-utils.ts went unnoticed. CI reported neither.
  eslint: {
    ignoreDuringBuilds: false,
  },
  typescript: {
    ignoreBuildErrors: false,
  },
  images: {
    unoptimized: true
  },
  // Ensure static assets are properly served
  assetPrefix: '',
  distDir: '.next',
  webpack: (config, { isServer }) => {
    if (!isServer) {
      config.resolve.fallback = {
        ...config.resolve.fallback,
        fs: false,
        net: false,
        tls: false,
      }
    }
    return config
  }
}

export default nextConfig
